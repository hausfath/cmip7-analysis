"""Weekly refresh, stage 3: mechanical checks on a regenerated analysis.

Compares the regenerated outputs with those in the last commit (git HEAD) and lists
every reason the refresh needs a person to look at it before it is published. An
empty list means the new data fit rules already documented in METHODS.md.

Writes <logdir>/check.json ({"clean": bool, "reasons": [...], "notes": [...], "diffs": {...}}).
Usage: python -I scripts/weekly/check.py <logdir>
"""
import io, json, pathlib, subprocess, sys
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs"
# Models with a documented exception to the cos(lat) weighting check (METHODS §1.1)
COSLAT_EXCEPTIONS = {"FGOALS-g3-5"}
COSLAT_TOL = 0.02          # K or W m-2; all other models are <= 0.011
OFFSET_TOL = 0.5           # degC from the median offset_5yr over all joins
# Changes to existing rows large enough to need a look
TOL = {"ECS_rm21": 0.2, "dT_2012_2021": 0.1, "dT_2081_2100": 0.3}
# Join rules other than "metadata" that were accepted for one model (METHODS §1.1); new runs of
# that model under the same rule are reported as notes, any other model needs review
ACCEPTED_RULES = {("FGOALS-g3-5", "relabel-same-label"), ("KACE-2-0-G", "placeholder-branch")}
# Members whose scenario joins are kept despite an unusual offset (user decision 2026-11-24; METHODS §1.1)
ACCEPTED_OFFSETS = {("EC-Earth3-ESM-1-1", "r4i1p1f1")}


def head_csv(rel):
    r = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=ROOT, capture_output=True, text=True)
    return pd.read_csv(io.StringIO(r.stdout)) if r.returncode == 0 else None


def main(logdir):
    logdir = pathlib.Path(logdir)
    scan = json.loads((logdir / "scan.json").read_text())
    reasons, notes, diffs = [], [], {}

    # 1. ESGF-side changes that are not simple additions
    if scan["new_models"]:
        reasons.append(f"New model(s): {', '.join(scan['new_models'])}. Needs holdings rows, model notes "
                       "and a check of figure colours and report text.")
    for k, label in [("datasets_removed", "dataset(s) withdrawn from ESGF"),
                     ("datasets_reversioned", "dataset(s) with a new version (refetched)"),
                     ("files_removed", "file(s) withdrawn from ESGF"),
                     ("checksum_changed", "file(s) with a changed checksum (refetched)")]:
        if scan[k]:
            reasons.append(f"{len(scan[k])} {label}: " + "; ".join(map(str, scan[k][:5]))
                           + (" ..." if len(scan[k]) > 5 else ""))

    # 2. Acquisition complete, with no new stale-checksum or unapproved TLS exceptions
    m = pd.read_csv(ROOT / "inventory" / "manifest_latest.csv")
    done = [(ROOT / "derived" / "monthly" / k / (f[:-3] + ".csv")).exists() for k, f in zip(m.dskey, m.filename)]
    missing = m.loc[[not d for d in done]]
    if len(missing):
        reasons.append(f"{len(missing)} manifest file(s) not acquired ("
                       + ", ".join(sorted(missing.source_id.unique())) + "); see acquire log.")
    log = pd.read_csv(ROOT / "derived" / "acquire_log.csv")
    log0 = head_csv("derived/acquire_log.csv")
    newlog = log.iloc[len(log0):] if log0 is not None else log
    if (newlog.status == "ok-stale").any():
        reasons.append(f"{(newlog.status == 'ok-stale').sum()} file(s) accepted with a stale index checksum.")
    tls = newlog[newlog.msg.astype(str).str.contains("tls-unverified")]
    if len(tls):
        notes.append(f"{len(tls)} file(s) fetched with TLS checks off for esg-earthlab.iap.ac.cn "
                     "(exact sha256 match required, per the documented exception).")

    # 3. Area weighting check on newly acquired datasets
    for p in (ROOT / "derived" / "attrs").glob("*.json"):
        if subprocess.run(["git", "cat-file", "-e", f"HEAD:{p.relative_to(ROOT)}"], cwd=ROOT,
                          capture_output=True).returncode == 0:
            continue
        a = json.loads(p.read_text())
        src = p.name.split(".")[4]
        d = a.get("coslat_minus_bnds_first_step")
        if d is not None and abs(float(d)) > COSLAT_TOL and src not in COSLAT_EXCEPTIONS:
            reasons.append(f"cos(lat) weighting check {float(d):+.3f} for {p.stem} (tolerance {COSLAT_TOL}).")

    # 4. Scenario joins: new joins must follow the metadata rule and look ordinary
    sc, sc0 = pd.read_csv(OUT / "splice_check.csv"), head_csv("outputs/splice_check.csv")
    key = ["model", "scenario", "member"]
    j = sc.merge(sc0[key + ["rule"]], on=key, how="left", suffixes=("", "_old"), indicator=True)
    new_j = j[j._merge == "left_only"]
    med = sc.loc[sc.used, "offset_5yr"].median()
    for r in new_j.itertuples():
        if (r.model, r.rule) in ACCEPTED_RULES:
            notes.append(f"New scenario join {r.model} {r.scenario} {r.member} under the documented "
                         f"'{r.rule}' rule (offset_5yr {r.offset_5yr:+.2f} degC).")
        elif r.rule != "metadata":
            reasons.append(f"New scenario join {r.model} {r.scenario} {r.member}: rule '{r.rule}'"
                           f" (parent found: {r.parent_found}).")
        elif pd.notna(r.offset_5yr) and abs(r.offset_5yr - med) > OFFSET_TOL:
            (notes if (r.model, r.member) in ACCEPTED_OFFSETS else reasons).append(f"New scenario join {r.model} {r.scenario} {r.member}: offset_5yr "
                           f"{r.offset_5yr:+.2f} degC vs median {med:+.2f}.")
    old_j = j[j._merge == "both"]
    for r in old_j[old_j.rule != old_j.rule_old].itertuples():
        if r.rule_old == "rejected" and r.rule == "metadata":   # parent now published
            notes.append(f"Join {r.model} {r.scenario} {r.member} now used (parent published).")
            continue
        reasons.append(f"Existing join {r.model} {r.scenario} {r.member} changed rule {r.rule_old} -> {r.rule}.")
    diffs["new_joins"] = len(new_j)

    # 5. Metric tables: new models need review; large changes to existing rows need review
    for f, k, col in [("cmip7_sensitivity.csv", ["model", "member"], "ECS_rm21"),
                      ("historical_metrics.csv", ["model", "track"], "dT_2012_2021"),
                      ("projection_metrics.csv", ["model", "scenario"], "dT_2081_2100")]:
        a, b = pd.read_csv(OUT / f), head_csv(f"outputs/{f}")
        new_models = sorted(set(a.model) - set(b.model))
        if new_models and f != "projection_metrics.csv":
            reasons.append(f"{f}: model(s) appear for the first time: {', '.join(new_models)}.")
        mm = a.merge(b, on=k, suffixes=("", "_old")).copy()
        mm["delta"] = mm[col] - mm[f"{col}_old"]
        big = mm[mm.delta.abs() > TOL[col]]
        for r in big.itertuples():
            reasons.append(f"{f}: {' '.join(str(getattr(r, c)) for c in k)} {col} changed by {r.delta:+.2f}.")
        added = a.merge(b[k], on=k, how="left", indicator=True)
        diffs[f] = {"rows_added": [" ".join(map(str, t)) for t in added.loc[added._merge == "left_only", k].itertuples(index=False)],
                    "rows_changed": mm.loc[mm.delta.abs() > 1e-9, k + [col, "delta"]].round(4).to_dict("records")}

    # 6. Holdings (member counts of tas datasets per model x experiment) for the docs tables
    def holdings(man):
        t = man[man.variable_id == "tas"].drop_duplicates("dskey")
        return t.groupby(["source_id", "experiment_id"]).variant_label.nunique()
    m0 = pd.read_csv(ROOT / scan["baseline_manifest"])
    h1, h0 = holdings(m), holdings(m0)
    hd = pd.concat([h0.rename("before"), h1.rename("after")], axis=1).fillna(0).astype(int)
    diffs["holdings_changed"] = hd[hd.before != hd.after].reset_index().to_dict("records")
    diffs["totals"] = {"files": len(m), "datasets": int(m.dskey.nunique()),
                       "GB": round(m["size"].sum() / 1e9, 1), "models": int(m.source_id.nunique())}

    res = {"clean": not reasons, "reasons": reasons, "notes": notes, "diffs": diffs}
    (logdir / "check.json").write_text(json.dumps(res, indent=1, default=str))
    print("CLEAN" if res["clean"] else "NEEDS REVIEW")
    for x in reasons:
        print(" -", x)
    for x in notes:
        print(" note:", x)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
