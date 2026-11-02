"""Stream CMIP7 monthly files from ESGF, reduce each to an area-weighted global
mean time series, and delete the raw file.

Never holds more than one raw file per worker on disk. Resumable: a file whose
reduced CSV already exists is skipped.

Inputs : inventory/manifest_latest.csv            (from build_manifest.py)
Outputs: derived/monthly/<dskey>/<filename>.csv   (one per raw file)
         derived/attrs/<dskey>.json               (global attrs of first file)
         derived/acquire_log.csv                  (status per file)

Usage  : python -I scripts/acquire_global_means.py [--sources A,B] [--exclude C]
"""
import argparse, hashlib, json, pathlib, sys, threading, time, traceback, urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests
import netCDF4
import cftime

ROOT = pathlib.Path(__file__).resolve().parents[1]
MAN = ROOT / "inventory" / "manifest_latest.csv"
OUT = ROOT / "derived" / "monthly"
ATTR = ROOT / "derived" / "attrs"
RAW = ROOT / "raw_tmp"          # transient; each file deleted after reduction
LOG = ROOT / "derived" / "acquire_log.csv"

HOST_LIMIT = {"esgf.ceda.ac.uk": 6, "crd-esgf-drc.ec.gc.ca": 5, "esgf.nci.org.au": 4}
DEFAULT_HOST_LIMIT = 2
ATTR_KEYS = ["branch_time_in_child", "branch_time_in_parent", "parent_experiment_id",
             "parent_variant_label", "parent_source_id", "parent_time_units", "parent_mip_era",
             "mip_era", "creation_date", "source", "grid", "nominal_resolution", "license",
             "license_id", "institution", "tracking_id", "variant_label", "version",
             "forcing_index", "physics_index", "external_variables", "title"]

_sems = {}
_log_lock = threading.Lock()


def host_sem(url):
    h = urllib.parse.urlparse(url).netloc
    if h not in _sems:
        _sems[h] = threading.Semaphore(HOST_LIMIT.get(h, DEFAULT_HOST_LIMIT))
    return _sems[h]


def verify(path, checksum, ctype):
    """ESGF CMIP7 checksums are sha256 multihashes ('1220' + hex digest)."""
    if not isinstance(checksum, str) or not checksum:
        return "unchecked"
    want = checksum[4:] if checksum.startswith("1220") and len(checksum) == 68 else checksum
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for ch in iter(lambda: f.read(1 << 22), b""):
            h.update(ch)
    if h.hexdigest() != want.lower():
        raise IOError(f"checksum mismatch {path.name}")
    return "sha256-ok"


ACCEPT_STALE = False   # --accept-stale: see process()


class StaleIndex(Exception):
    """Server file is complete (matches its own Content-Length) but differs from the index."""


def download(urls, dest, size):
    last = None
    for attempt in range(4):
        for url in urls:
            try:
                with host_sem(url):
                    with requests.get(url, stream=True, timeout=(30, 300)) as r:
                        r.raise_for_status()
                        clen = int(r.headers.get("content-length", -1))
                        with open(dest, "wb") as f:
                            for ch in r.iter_content(1 << 22):
                                f.write(ch)
                got = dest.stat().st_size
                if size and got != int(size):
                    if ACCEPT_STALE and got == clen:
                        raise StaleIndex(url)
                    raise IOError(f"size {got} != {size}")
                return url
            except StaleIndex:
                raise
            except Exception as e:  # noqa: BLE001
                last = e
                dest.unlink(missing_ok=True)
                time.sleep(5 * (attempt + 1))
    raise last


def area_weights(ds):
    """Weights for a regular lon-lat grid from latitude bounds (exact spherical
    band areas); falls back to cos(lat) only if no bounds are present."""
    lat = ds.variables["lat"]
    if lat.ndim != 1:
        raise ValueError("non-rectilinear grid: needs areacella")
    bname = getattr(lat, "bounds", None)
    if bname and bname in ds.variables:
        b = np.asarray(ds.variables[bname][:], dtype=float)
        w = np.abs(np.sin(np.deg2rad(b[:, 1])) - np.sin(np.deg2rad(b[:, 0])))
        kind = "lat_bnds"
    else:
        w = np.cos(np.deg2rad(np.asarray(lat[:], dtype=float)))
        kind = "coslat"
    return w, kind


def reduce_file(path, var):
    with netCDF4.Dataset(path) as ds:
        attrs = {k: str(ds.getncattr(k)) for k in ATTR_KEYS if k in ds.ncattrs()}
        t = ds.variables["time"]
        units, cal = t.units, getattr(t, "calendar", "standard")
        tv = np.asarray(t[:], dtype=float)
        dates = cftime.num2date(tv, units, cal)
        wlat, kind = area_weights(ds)
        v = ds.variables[var]
        dims = v.dimensions
        assert v.ndim == 3 and dims[0] == "time" and dims[-2:] == ("lat", "lon"), dims
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        assert np.allclose(np.diff(lon), np.diff(lon)[0]), "non-uniform lon spacing"
        out = np.empty(len(tv))
        step = 120
        for i in range(0, len(tv), step):
            x = np.ma.masked_invalid(np.ma.asarray(v[i:i + step], dtype=float))
            W = np.broadcast_to(wlat[None, :, None], x.shape)
            W = np.ma.array(W, mask=np.ma.getmaskarray(x))
            out[i:i + step] = (x * W).sum(axis=(1, 2)) / W.sum(axis=(1, 2))
            nmask = int(np.ma.getmaskarray(x).sum())
            if nmask:
                raise ValueError(f"{nmask} masked/NaN points in a global field")
        # cos(lat) alternative on the first time step, for the weighting check
        x0 = np.ma.masked_invalid(np.ma.asarray(v[0], dtype=float)).reshape(-1, len(wlat), v.shape[-1])[0]
        c = np.cos(np.deg2rad(np.asarray(ds.variables["lat"][:], dtype=float)))
        attrs["coslat_minus_bnds_first_step"] = str(float((x0 * c[:, None]).sum() / (c.sum() * x0.shape[1])) - out[0])
        attrs.update(time_units=units, calendar=cal, weight_kind=kind,
                     nlat=str(len(wlat)), nlon=str(v.shape[-1]), units=getattr(v, "units", ""))
    df = pd.DataFrame({"time": tv, "year": [d.year for d in dates], "month": [d.month for d in dates],
                       var: out})
    return df, attrs


def process(row):
    dskey, fn, var = row.dskey, row.filename, row.variable_id
    odir = OUT / dskey
    ocsv = odir / (fn[:-3] + ".csv")
    if ocsv.exists():
        return dskey, fn, "skip", ""
    odir.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(exist_ok=True)
    tmp = RAW / f"{dskey}__{fn}"
    t0 = time.time()
    try:
        try:
            url = download(json.loads(row.urls), tmp, row.size)
            chk = verify(tmp, row.checksum, row.checksum_type)
        except StaleIndex as e:
            url, chk = str(e), "index-stale: size==server Content-Length, sha256!=index"
        df, attrs = reduce_file(tmp, var)
        df.to_csv(ocsv.with_suffix(".part"), index=False)
        ocsv.with_suffix(".part").rename(ocsv)
        af = ATTR / f"{dskey}.json"
        if not af.exists():
            ATTR.mkdir(parents=True, exist_ok=True)
            attrs.update(url=url, checksum_status=chk)
            af.write_text(json.dumps(attrs, indent=1))
        status = "ok-stale" if chk.startswith("index-stale") else "ok"
        msg = f"{row.size/1e6/(time.time()-t0):.1f} MB/s {chk}"
    except Exception as e:  # noqa: BLE001
        status, msg = "fail", f"{type(e).__name__}: {e}"
        traceback.print_exc()
    finally:
        tmp.unlink(missing_ok=True)
    with _log_lock:
        new = not LOG.exists()
        with open(LOG, "a") as f:
            if new:
                f.write("utc,dskey,filename,status,msg\n")
            f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())},{dskey},{fn},{status},\"{msg}\"\n")
    return dskey, fn, status, msg


# Order: what ECS/TCR need first, then historical, then scenarios; the
# large high-resolution model last.
PRIO = {"abrupt-4xCO2": 0, "piControl": 0, "1pctCO2": 0, "historical": 1, "esm-hist": 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", default="")
    ap.add_argument("--exclude", default="")
    ap.add_argument("--accept-stale", action="store_true",
                    help="accept files whose size matches the server Content-Length but not the index")
    a = ap.parse_args()
    global ACCEPT_STALE
    ACCEPT_STALE = a.accept_stale
    m = pd.read_csv(MAN)
    if a.sources:
        m = m[m.source_id.isin(a.sources.split(","))]
    if a.exclude:
        m = m[~m.source_id.isin(a.exclude.split(","))]
    m = m.assign(p=m.experiment_id.map(PRIO).fillna(2) + (m.source_id == "UKCM2a-0-HH") * 3)
    # One pool per data node, sized to that node's connection cap, so a slow
    # node never idles workers that could be pulling from another node.
    m["host"] = m.urls.map(lambda u: urllib.parse.urlparse(json.loads(u)[0]).netloc)
    m = m.sort_values(["p", "dskey", "filename"])
    todo = [r for r in m.itertuples() if not (OUT / r.dskey / (r.filename[:-3] + ".csv")).exists()]
    print(f"{len(todo)} files to fetch, {sum(r.size for r in todo)/1e9:.1f} GB", flush=True)
    done = gb = 0
    t0 = time.time()
    pools, futs = [], {}
    for host in sorted({r.host for r in todo}):
        ex = ThreadPoolExecutor(HOST_LIMIT.get(host, DEFAULT_HOST_LIMIT))
        pools.append(ex)
        futs.update({ex.submit(process, r): r for r in todo if r.host == host})
    for f in as_completed(futs):
        r = futs[f]
        dskey, fn, status, msg = f.result()
        done += 1
        gb += r.size / 1e9 if status == "ok" else 0
        if done % 25 == 0 or status == "fail":
            el = (time.time() - t0) / 3600
            print(f"[{done}/{len(todo)}] {gb:.1f} GB in {el:.2f} h ({gb/el if el else 0:.1f} GB/h) "
                  f"last: {status} {fn} {msg}", flush=True)
    for ex in pools:
        ex.shutdown()
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
