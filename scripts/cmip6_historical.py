"""CMIP6 historical GSAT metrics from the Hausfather et al. (2022) supporting data (METHODS.md §2.9).

  python -I scripts/cmip6_historical.py

One annual series per model (single member, inferred from interannual noise; the member and
variable are not documented in H22 or Zenodo 6476375). Historical (1850–2014) is spliced to
SSP2-4.5 (2015–) as the primary; SSP1-2.6 / SSP3-7.0 / SSP5-8.5 are the splice sensitivity.
Series are re-baselined from H22's 1850–1899 to 1850–1900 to match the obs and CMIP7 rows.

Outputs
  outputs/cmip6_historical_metrics.csv     per model (primary splice), plus the obs row from
                                           historical_metrics.csv
  outputs/cmip6_splice_sensitivity.csv     per model × SSP metrics and 2014→2015 jump
  outputs/cmip6_h22_checks.txt             read check, baseline shift, splice checks, medians
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from gmst import METRICS, metrics, BASE  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
XL = ROOT / "reference" / "hausfather2022_supporting_data.xlsx"
OUT = ROOT / "outputs"
SSPS = {"SSP245": "SSP2-4.5", "SSP126": "SSP1-2.6", "SSP370": "SSP3-7.0", "SSP585": "SSP5-8.5"}
PRIMARY = "SSP245"
# CMIP6 predecessors of the CMIP7 models (figures_sensitivity.LINEAGE); HadGEM3-GC31-MM is not in H22
PREDECESSOR = {"CanESM5": "CanESM6-0-MR, CanESM5-1", "UKESM1-0-LL": "UKESM1-3-LL",
               "HadGEM3-GC31-LL": "UKCM2-0-LL", "ACCESS-ESM1-5": "ACCESS-ESM1-6"}


def sheet(name):
    d = pd.read_excel(XL, name, index_col=0)
    d.index = d.index.astype(int)
    return d


def main():
    hist = sheet("Historical")
    ssp = {k: sheet(k) for k in SSPS}
    log = []

    # Read check: the workbook's "CMIP6 MMM" central column (labelled 50th percentile) should be the
    # across-model mean of the SSP sheet columns. Its 5th/95th columns match no simple statistic of
    # the 37 columns (checked: quantiles, min/max, mean ± 1.645 sd), so they are not used.
    mmm = sheet("CMIP6 MMM")
    for k, lab in SSPS.items():
        d = (ssp[k].mean(axis=1) - mmm[f"{lab} 50th Percentile"]).abs()
        log.append(f"{lab}: max |37-model mean − MMM '50th Percentile'| 2015–2099 = {d.loc[:2099].max():.1e}")
    # Single member vs ensemble mean: interannual noise (SD of 1850–1950 first differences / √2)
    noise = hist.loc[1850:1950].diff().std() / np.sqrt(2)
    log.append(f"Interannual noise 1850–1950: {noise.min():.3f}–{noise.max():.3f} °C across models "
               f"(CanESM5 {noise['CanESM5']:.3f}, MIROC6 {noise['MIROC6']:.3f}; a ~50-member mean would be "
               f"~{noise.median() / np.sqrt(50):.3f})")
    log.append(f"H22 1850–1899 mean, max |·| over models = {hist.loc[1850:1899].mean().abs().max():.1e}")
    shift = hist.loc[BASE[0]:BASE[1]].mean()
    log.append(f"Re-baselining 1850–1899 → 1850–1900 shifts models by {shift.min():+.4f} to {shift.max():+.4f} °C")

    rows, sens = [], []
    for m in hist.columns:
        for k in SSPS:
            s = ssp[k][m] if m in ssp[k] else None
            if s is None or s.loc[2015:2024].isna().any():
                sens.append({"model": m, "ssp": SSPS[k], "available": False})
                continue
            x = pd.concat([hist[m], s]) - shift[m]
            met = metrics(x)
            # 2014→2015 join vs the model's own 1985–2014 year-to-year changes
            dd = hist[m].loc[1985:2014].diff().dropna()
            jump = s.loc[2015] - hist[m].loc[2014]
            sens.append({"model": m, "ssp": SSPS[k], "available": True, "jump_2014_2015": jump,
                         "jump_z": (jump - dd.mean()) / dd.std(), **met})
            if k == PRIMARY:
                rows.append({"model": m, "track": "historical", "n_members": 1, "splice": SSPS[k],
                             "predecessor_of": PREDECESSOR.get(m, ""), **met})
    sens = pd.DataFrame(sens)
    hm = pd.DataFrame(rows).set_index("model")
    # Bar in the figure: range across the available SSP splices (stands in for fig02's member range)
    ok = sens[sens.available]
    for k in METRICS:
        g = ok.groupby("model")[k]
        hm[f"{k}_min"], hm[f"{k}_max"] = g.min(), g.max()
    hm["n_splices"] = ok.groupby("model").size()
    hm = hm.reset_index()

    # Splice sensitivity: spread across SSPs per model, for metrics ending by 2021 and by 2024
    for k in METRICS:
        sp = ok.groupby("model")[k].agg(lambda v: v.max() - v.min())
        log.append(f"splice spread across SSPs, {k}: max {sp.max():.4f} ({sp.idxmax()}), median {sp.median():.4f}")
    j = ok[ok.ssp == SSPS[PRIMARY]]
    log.append(f"2014→2015 jump ({SSPS[PRIMARY]}): median {j.jump_2014_2015.median():+.3f}, range "
               f"{j.jump_2014_2015.min():+.3f} to {j.jump_2014_2015.max():+.3f} °C; |z| > 2 for "
               f"{(j.jump_z.abs() > 2).sum()} of {len(j)} (≈{0.0455 * len(j):.1f} expected by chance)")
    log.append("Missing SSPs: " + "; ".join(f"{r.model} {r.ssp}" for r in sens[~sens.available].itertuples()))

    # Sanity: across-model median warming vs obs
    obs = pd.read_csv(OUT / "historical_metrics.csv")
    obs = obs[obs.track == "obs"]
    for k in ["dT_2011_2020", "dT_2012_2021", "trend_1980_2021"]:
        log.append(f"{k}: CMIP6 median {hm[k].median():.3f}, 5–95% {hm[k].quantile(.05):.3f}–"
                   f"{hm[k].quantile(.95):.3f}; obs median {obs[k].iloc[0]:.3f}")

    pd.concat([obs, hm], ignore_index=True).to_csv(OUT / "cmip6_historical_metrics.csv", index=False)
    sens.to_csv(OUT / "cmip6_splice_sensitivity.csv", index=False)
    (OUT / "cmip6_h22_checks.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
