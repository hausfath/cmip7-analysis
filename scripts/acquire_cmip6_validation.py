"""Acquire CMIP6 global-mean monthly tas/rsdt/rsut/rlut for a few models from the
Pangeo cloud archive (gs://cmip6) so the Gregory ECS / TCR code can be validated
against published values (Zelinka et al. 2020 v2.2 / Hausfather et al. 2022 SI)
before it is applied to CMIP7.

Run in a separate venv with zarr<3 and gcsfs (requirements-acq.txt):
  .venv-acq/bin/python scripts/acquire_cmip6_validation.py

Output: derived/cmip6_validation/<model>.<exp>.<member>.<var>.csv  (time,year,month,<var>)
        derived/cmip6_validation/<model>.<exp>.<member>.attrs.json
"""
import json, pathlib, sys, warnings
import numpy as np
import pandas as pd
import xarray as xr

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "derived" / "cmip6_validation"
OUT.mkdir(parents=True, exist_ok=True)
CAT = "https://storage.googleapis.com/cmip6/pangeo-cmip6.csv"

# (model, abrupt/1pct member, piControl member). The first three are the code
# validation set; the rest are CMIP6 lineage predecessors of CMIP7 models, run
# through the same code so CMIP6->CMIP7 differences are method-consistent.
RUNS = [("CanESM5", "r1i1p1f1", "r1i1p1f1"), ("UKESM1-0-LL", "r1i1p1f2", "r1i1p1f2"),
        ("ACCESS-ESM1-5", "r1i1p1f1", "r1i1p1f1"), ("CanESM5", "r1i1p2f1", "r1i1p2f1"),
        ("HadGEM3-GC31-LL", "r1i1p1f3", "r1i1p1f1"), ("HadGEM3-GC31-MM", "r1i1p1f3", "r1i1p1f1"),
        ("CNRM-ESM2-1", "r1i1p1f2", "r1i1p1f2")]
EXPS = {"abrupt-4xCO2": ["tas", "rsdt", "rsut", "rlut"], "piControl": ["tas", "rsdt", "rsut", "rlut"],
        "1pctCO2": ["tas"]}
KEEP = ["branch_time_in_parent", "branch_time_in_child", "parent_variant_label", "parent_time_units",
        "parent_experiment_id", "version_id", "tracking_id"]


def gmean(da):
    lat = da["lat"]
    # Pangeo zarr usually drops bounds; use exact band areas from midpoints
    l = lat.values.astype(float)
    edges = np.concatenate([[-90.0], 0.5 * (l[1:] + l[:-1]), [90.0]])
    w = xr.DataArray(np.abs(np.sin(np.deg2rad(edges[1:])) - np.sin(np.deg2rad(edges[:-1]))),
                     dims="lat", coords={"lat": lat})
    return da.weighted(w).mean(("lat", "lon"))


def write_attrs(model, exp, mem, var, ds, url):
    """Merge into the per-run attrs file. Time units are stored per variable
    (time_units_<var>): stores of the same run can differ (CanESM5 r1i1p1f1
    abrupt-4xCO2 tas uses 'hours since 1850-01-16 12:00', the other variables
    'days since 1850-01-01')."""
    f = OUT / f"{model}.{exp}.{mem}.attrs.json"
    a = json.loads(f.read_text()) if f.exists() else {}
    a.update({k: str(ds.attrs.get(k)) for k in KEEP})
    tu, cal = ds["time"].attrs["units"], ds["time"].attrs.get("calendar", "standard")
    a.update({"time_units": a.get("time_units", tu), f"time_units_{var}": tu, "calendar": cal,
              f"zstore_{var}": url})
    f.write_text(json.dumps(a, indent=1))


def backfill(cat):
    """Metadata-only pass recording per-variable time units for existing CSVs."""
    for model, cmem, pmem in RUNS:
        for exp, vars_ in EXPS.items():
            mem = pmem if exp == "piControl" else cmem
            for var in vars_:
                if not (OUT / f"{model}.{exp}.{mem}.{var}.csv").exists():
                    continue
                q = cat[(cat.source_id == model) & (cat.experiment_id == exp) & (cat.member_id == mem)
                        & (cat.table_id == "Amon") & (cat.variable_id == var)]
                url = q.sort_values("version").zstore.iloc[-1]
                ds = xr.open_zarr(url, consolidated=True, storage_options={"token": "anon"}, decode_times=False)
                write_attrs(model, exp, mem, var, ds, url)
                print("attrs", model, exp, var, ds["time"].attrs["units"], flush=True)


def fetch_one(cat, model, exp, mem, var):
    out = OUT / f"{model}.{exp}.{mem}.{var}.csv"
    q = cat[(cat.source_id == model) & (cat.experiment_id == exp) & (cat.member_id == mem)
            & (cat.table_id == "Amon") & (cat.variable_id == var)]
    if q.empty:
        print("MISSING", model, exp, mem, var, flush=True)
        return
    url = q.sort_values("version").zstore.iloc[-1]
    ds = xr.open_zarr(url, consolidated=True, storage_options={"token": "anon"}, decode_times=False)
    g = gmean(ds[var]).load()
    tu, cal = ds["time"].attrs["units"], ds["time"].attrs.get("calendar", "standard")
    import cftime
    dates = cftime.num2date(g["time"].values, tu, cal)
    tmp = out.with_suffix(".part")
    pd.DataFrame({"time": g["time"].values, "year": [d.year for d in dates],
                  "month": [d.month for d in dates], var: g.values}).to_csv(tmp, index=False)
    tmp.rename(out)
    write_attrs(model, exp, mem, var, ds, url)
    print("ok", model, exp, var, len(g), flush=True)


def main():
    cat = pd.read_csv(CAT)
    if "--backfill-attrs" in sys.argv:
        return backfill(cat)
    if "--one" in sys.argv:
        return fetch_one(cat, *sys.argv[sys.argv.index("--one") + 1:][:4])
    # Each variable in a child process under a timeout: gcsfs reads can stall
    # indefinitely after a network drop (seen twice), and a stall must not block the run.
    import subprocess
    for model, cmem, pmem in RUNS:
        for exp, vars_ in EXPS.items():
            mem = pmem if exp == "piControl" else cmem
            for var in vars_:
                if (OUT / f"{model}.{exp}.{mem}.{var}.csv").exists():
                    continue
                for attempt in range(3):
                    try:
                        subprocess.run([sys.executable, __file__, "--one", model, exp, mem, var],
                                       timeout=900, check=True)
                        break
                    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
                        print("retry", model, exp, var, type(e).__name__, flush=True)


if __name__ == "__main__":
    main()
