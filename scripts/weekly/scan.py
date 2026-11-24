"""Weekly refresh, stage 1: compare a fresh ESGF manifest with the last committed one.

Run after scripts/build_manifest.py has written inventory/manifest_<today>.csv.
The baseline is the newest dated manifest in the last commit (git HEAD).

Also lists models that have DECK idealised runs (abrupt-4xCO2 / 1pctCO2) on ESGF
but not yet the variables needed for ECS/TCR (tas + rsdt/rsut/rlut in abrupt-4xCO2
and piControl), so partly published models are visible before they become usable.

Writes <logdir>/scan.json. Exit code 0 = nothing new, 10 = new or changed files.
Usage: python -I scripts/weekly/scan.py <logdir>
"""
import json, pathlib, subprocess, sys
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
INV = ROOT / "inventory"
NEED_TAS = ("tas", "tavg-h2m-hxy-u")
NEED_RAD = {("rsdt", "tavg-u-hxy-u"), ("rsut", "tavg-u-hxy-u"), ("rlut", "tavg-u-hxy-u")}


def committed_manifests():
    out = subprocess.run(["git", "ls-tree", "--name-only", "HEAD", "inventory/"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout.split()
    return sorted(p for p in out if pathlib.Path(p).name.startswith("manifest_20"))


def deck_coverage():
    """Models with abrupt-4xCO2/1pctCO2 datasets of any variable vs the ECS/TCR variables."""
    from intake_esgf import ESGFCatalog
    cat = ESGFCatalog()
    cat.search(project="CMIP7", experiment_id=["abrupt-4xCO2", "1pctCO2", "piControl"], frequency="mon")
    df = cat.df
    rows = []
    for src, g in df.groupby("source_id"):
        have = set(zip(g.experiment_id, g.variable_id, g.variable_branding_suffix))
        def ok(exp, need):
            return all((exp, v, b) in have for v, b in need)
        idealised = sorted(set(g.experiment_id) & {"abrupt-4xCO2", "1pctCO2"})
        if not idealised:
            continue
        rows.append({"model": src, "idealised_experiments": idealised,
                     "ecs_ready": ok("abrupt-4xCO2", NEED_RAD | {NEED_TAS}) and ok("piControl", NEED_RAD | {NEED_TAS}),
                     "tcr_ready": ok("1pctCO2", {NEED_TAS}) and ok("piControl", {NEED_TAS}),
                     "n_variables_abrupt": int(g[g.experiment_id == "abrupt-4xCO2"].variable_id.nunique())})
    return rows


def main(logdir):
    logdir = pathlib.Path(logdir)
    logdir.mkdir(parents=True, exist_ok=True)
    dated = sorted(INV.glob("manifest_20*.csv"))
    new_path = dated[-1]
    base_rel = [p for p in committed_manifests() if pathlib.Path(p).name != new_path.name][-1]
    new = pd.read_csv(new_path)
    base = pd.read_csv(ROOT / base_rel)

    ds_new = new.groupby("dskey").agg(version=("version", "first"), source_id=("source_id", "first"),
                                      experiment_id=("experiment_id", "first"), variable_id=("variable_id", "first"))
    ds_base = base.groupby("dskey").agg(version=("version", "first"))
    added = ds_new.index.difference(ds_base.index)
    removed = ds_base.index.difference(ds_new.index)
    common = ds_new.index.intersection(ds_base.index)
    reversioned = [k for k in common if ds_new.version[k] != ds_base.version[k]]
    fk_new = set(zip(new.dskey, new.filename))
    fk_base = set(zip(base.dskey, base.filename))
    files_added = sorted(fk_new - fk_base)
    files_removed = sorted(fk_base - fk_new)
    cn = new.set_index(["dskey", "filename"]).checksum
    cb = base.set_index(["dskey", "filename"]).checksum
    both = cn.index.intersection(cb.index)
    checksum_changed = [list(k) for k in both if cn[k] != cb[k]]

    summary = (ds_new.loc[added].groupby(["source_id", "experiment_id"]).size()
               .rename("datasets").reset_index().to_dict("records"))
    try:
        deck = deck_coverage()
    except Exception as e:  # the coverage listing is informational; never block the scan on it
        deck = {"error": repr(e)}
    res = {
        "baseline_manifest": base_rel, "new_manifest": str(new_path.relative_to(ROOT)),
        "baseline_files": len(base), "new_files": len(new),
        "new_GB": round(new["size"].sum() / 1e9, 1),
        "models_baseline": sorted(base.source_id.unique()), "models_new": sorted(new.source_id.unique()),
        "new_models": sorted(set(new.source_id) - set(base.source_id)),
        "datasets_added": len(added), "datasets_added_by_model_experiment": summary,
        "datasets_removed": sorted(removed), "datasets_reversioned": reversioned,
        "files_added": len(files_added), "files_added_GB": round(
            new.set_index(["dskey", "filename"]).loc[files_added, "size"].sum() / 1e9, 2) if files_added else 0.0,
        "files_removed": [list(k) for k in files_removed], "checksum_changed": checksum_changed,
        "deck_coverage": deck,
    }
    res["changed"] = bool(len(added) or len(removed) or reversioned or files_added or files_removed or checksum_changed)
    (logdir / "scan.json").write_text(json.dumps(res, indent=1, default=str))
    # acquire_global_means.py skips files whose reduced CSV exists, so files in a new
    # dataset version or with a new checksum must have their old CSV removed to be refetched
    stale = {(k, f) for k, f in fk_new if k in set(reversioned)} | {tuple(k) for k in checksum_changed}
    (logdir / "refetch.txt").write_text("".join(f"derived/monthly/{k}/{f[:-3]}.csv\n" for k, f in sorted(stale)))
    print(f"baseline {base_rel}: {len(base)} files; new: {len(new)} files; +{len(added)} datasets "
          f"(+{len(files_added)} files, {res['files_added_GB']} GB); removed {len(removed)}; "
          f"reversioned {len(reversioned)}; checksum changed {len(checksum_changed)}; new models {res['new_models']}")
    if isinstance(deck, list):
        for d in deck:
            if not d["ecs_ready"]:
                print(f"DECK not ECS-ready: {d['model']} {d['idealised_experiments']} "
                      f"({d['n_variables_abrupt']} variables in abrupt-4xCO2)")
    return 10 if res["changed"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
