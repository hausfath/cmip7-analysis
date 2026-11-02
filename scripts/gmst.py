"""CMIP7 historical and projected GSAT vs the GSAT ensemble (METHODS.md §2.6–2.8).

  python -I scripts/gmst.py                       obs from the vendored summaries in data/obs/
  python -I scripts/gmst.py --obs <ensemble.csv>  rebuild data/obs/ from the full 10,000-member
                                                  GSAT ensemble (not public), then run

Only the obs statistics this script uses are vendored: the annual median / 5 / 95 % / mean, and
the median, 5 / 95 % and SD of each metric over the members. Both modes give identical outputs.

Scenario runs are joined to the historical run named in their own parent metadata
(parent_experiment_id / parent_variant_label, branch at 2022-01-01), never by
matching variant labels: some UKESM1-3-LL esm-scen7 f1 members branch from f2
esm-hist members.

Outputs
  derived/gmst/<model>.<experiment>.csv   annual ΔT vs 1850–1900 (the parent historical
                                          run's own 1850–1900 mean), one column per member.
                                          Historical tracks are spliced to a scenario
                                          (preference m > h > vl) when one exists.
  outputs/obs_gsat_ensemble_summary.csv   obs median / 5 / 95 % per year
  outputs/historical_metrics.csv          per model × track, plus obs row
  outputs/projection_metrics.csv          per model × scenario
  outputs/splice_sensitivity.csv          ensemble-mean 2015–2024 warming and 1980–2024 trend with each
                                          scenario, over members continued in all of them
  outputs/splice_check.csv                parent lookup result for every scenario run
  outputs/drift_check.csv                 parent-control drift over each historical run
"""
import numpy as np
import pandas as pd
import cftime
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # -I drops the script dir
import cmiplib as L  # noqa: E402

ROOT = L.ROOT
OBS_DIR = ROOT / "data" / "obs"
GD = ROOT / "derived" / "gmst"
OUT = ROOT / "outputs"
BASE = (1850, 1900)
HIST_END = 2021
MODELS = ["CanESM6-0-MR", "CanESM5-1", "UKESM1-3-LL", "UKCM2-0-LL", "UKCM2a-0-HH", "ACCESS-ESM1-6",
          "CAS-ESM2-1-CSM", "EC-Earth3-ESM-1-1", "CNRM-ESM2-1e"]
SCEN = {"historical": ["scen7-m", "scen7-h", "scen7-vl"],              # splice preference order
        "esm-hist": ["esm-scen7-m", "esm-scen7-h", "esm-scen7-vl"]}
CONTROL = {"historical": "piControl", "esm-hist": "esm-piControl"}
# metric -> (kind, y0, y1); ΔT windows are relative to 1850–1900
METRICS = {"dT_2012_2021": ("mean", 2012, 2021), "dT_2015_2024": ("mean", 2015, 2024),
           "dT_2011_2020": ("mean", 2011, 2020),
           "trend_1980_2021": ("trend", 1980, 2021), "trend_1950_2021": ("trend", 1950, 2021),
           "trend_1980_2024": ("trend", 1980, 2024)}
AR6_2011_2020 = (1.09, 0.95, 1.20)   # AR6 WGI SPM A.1.2, GMST 2011–2020 vs 1850–1900


def metric(s, kind, y0, y1):
    s = s.loc[y0:y1].dropna()
    if len(s) < (y1 - y0 + 1):
        return np.nan
    if kind == "mean":
        return s.mean()
    return 10 * np.polyfit(s.index.to_numpy(float), s.to_numpy(), 1)[0]


def metrics(s):
    return {k: metric(s, *v) for k, v in METRICS.items()}


def obs_metrics(obs):
    """Vectorised over the 10,000 obs members."""
    out = {}
    for k, (kind, y0, y1) in METRICS.items():
        x = obs.loc[y0:y1]
        if len(x) < y1 - y0 + 1:
            out[k] = np.full(obs.shape[1], np.nan)
        elif kind == "mean":
            out[k] = x.mean().to_numpy()
        else:
            t = x.index.to_numpy(float) - x.index.to_numpy(float).mean()
            out[k] = 10 * (t[:, None] * x.to_numpy()).sum(0) / (t ** 2).sum()
    return pd.DataFrame(out)


def safe(source, exp, mem):
    try:
        return L.load_cmip7(source, exp, mem, "tas")
    except FileNotFoundError:
        return None, None


def find_scenarios(model, hexp):
    """{parent historical member: {scenario: (series, attrs)}} from scenario parent metadata."""
    out, checks = {}, []
    for s in SCEN[hexp]:
        for mem in L.members_cmip7(model, s):
            x, a = safe(model, s, mem)
            pexp, pmem = a.get("parent_experiment_id"), a.get("parent_variant_label")
            d = cftime.num2date(float(a["branch_time_in_parent"]), a.get("parent_time_units") or a["time_units"],
                                a["calendar"])
            ph, _ = safe(model, pexp, pmem)
            ok = pexp == hexp and ph is not None and (d.year, d.month, d.day) == (HIST_END + 1, 1, 1)
            jump = float(x.loc[HIST_END + 1] - ph.loc[HIST_END]) if ph is not None else np.nan
            checks.append({"model": model, "scenario": s, "member": mem, "parent_experiment": pexp,
                           "parent_member": pmem, "label_differs": pmem != mem,
                           "parent_found": ph is not None, "branch_decoded": str(d),
                           "jump_2021_2022": jump, "used": ok})
            if ok:
                out.setdefault(pmem, {})[s] = x
    return out, checks


def drift(model, hexp, mem, hattrs, pi_cache):
    """Parent-control OLS slope over the run's contemporaneous 1850–2021 segment, and the
    implied spurious contribution to the 2012–2021 minus 1850–1900 difference."""
    cexp = CONTROL[hexp]
    if cexp not in pi_cache:
        pmems = L.members_cmip7(model, cexp)
        pi_cache[cexp] = safe(model, cexp, pmems[0]) if pmems else (None, None)
    pi, pa = pi_cache[cexp]
    r = {"model": model, "track": hexp, "member": mem, "control": cexp}
    if pi is None:
        return {**r, "note": "no control run published"}
    by, off, _ = L.branch_year(hattrs, pa)
    n = HIST_END - BASE[0] + 1
    seg = pi.loc[by:by + n - 1]
    full = 10 * np.polyfit(pi.index.to_numpy(float), pi.to_numpy(), 1)[0]
    r.update(branch_year=by, branch_offset_days=off, control_years=f"{int(pi.index.min())}-{int(pi.index.max())}",
             control_trend_full_per_dec=full)
    if len(seg) == n:
        sl = 10 * np.polyfit(seg.index.to_numpy(float), seg.to_numpy(), 1)[0]
        r.update(control_trend_segment_per_dec=sl, segment="contemporaneous")
    else:
        sl = full
        r.update(control_trend_segment_per_dec=np.nan,
                 segment=f"only {len(seg)}/{n} contemporaneous years; full-record slope used")
    # The segment slope mixes drift with centennial variability; the full-record slope
    # isolates the long-term drift. Both reported.
    lever = np.mean([2012, 2021]) - np.mean(BASE)
    r["drift_in_dT_2012_2021"] = sl / 10 * lever
    r["drift_in_dT_2012_2021_fullrecord"] = full / 10 * lever
    return r


def build_model(model, pi_cache):
    res, splice_rows, checks, drifts, ivs = {}, [], [], [], {}
    for hexp in SCEN:
        mems = L.members_cmip7(model, hexp)
        if not mems:
            continue
        scen_by_parent, ch = find_scenarios(model, hexp)
        checks += ch
        hist_cols, used = {}, {}
        scen_cols = {s: {} for s in SCEN[hexp]}
        for mem in mems:
            h, ha = safe(model, hexp, mem)
            if h is None or h.index.min() > BASE[0]:
                continue
            h = h.loc[:HIST_END]
            base = h.loc[BASE[0]:BASE[1]].mean()
            drifts.append(drift(model, hexp, mem, ha, pi_cache))
            ext = {s: x.loc[HIST_END + 1:] for s, x in scen_by_parent.get(mem, {}).items()}
            for s, x in ext.items():
                scen_cols[s][mem] = pd.concat([h, x]) - base
            pick = next((s for s in SCEN[hexp] if s in ext), None)
            hist_cols[mem] = (pd.concat([h, ext[pick]]) if pick else h) - base
            used[mem] = pick or "none"
        # Splice sensitivity at the ensemble level, over members continued in every scenario
        # (per-member differences are dominated by chaotic divergence after 2022).
        have = [s for s in SCEN[hexp] if scen_cols[s]]
        common = sorted(set.intersection(*(set(scen_cols[s]) for s in have))) if len(have) > 1 else []
        if common:
            for s in have:
                mm = pd.DataFrame([metrics(scen_cols[s][m]) for m in common])
                splice_rows.append({"model": model, "track": hexp, "scenario": s, "n_common": len(common),
                                    "dT_2015_2024": mm.dT_2015_2024.mean(),
                                    "trend_1980_2024": mm.trend_1980_2024.mean()})
        if hist_cols:
            res[hexp] = (pd.DataFrame(hist_cols), used)
        for s, cols in scen_cols.items():
            if cols:
                res[s] = (pd.DataFrame(cols), None)
        pi = pi_cache.get(CONTROL[hexp], (None, None))[0]
        if pi is not None:
            ivs[hexp] = L.pi_variability(pi, lengths=(10, 51), trend_lengths=(42, 45, 72))
    return res, splice_rows, checks, drifts, ivs


def iv_sigma(metric_name, iv):
    """piControl σ for a metric (mean-window differences treat the 10-yr and 51-yr means
    as independent, ~140 yr apart)."""
    if iv is None:
        return np.nan
    kind, y0, y1 = METRICS[metric_name]
    if kind == "mean":
        return np.hypot(iv["sd_mean10"], iv["sd_mean51"])
    return iv[f"sd_trend{y1 - y0 + 1}"]


def vendor_obs(path):
    """Summarise the full obs ensemble (year x member, no header) into data/obs/."""
    obs = pd.read_csv(path, header=None, index_col=0)
    obs.index.name = "year"
    OBS_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"median": obs.median(axis=1), "p05": obs.quantile(0.05, axis=1),
                  "p95": obs.quantile(0.95, axis=1), "mean": obs.mean(axis=1)}).to_csv(OBS_DIR / "annual_summary.csv")
    om = obs_metrics(obs)
    st = pd.DataFrame({"median": om.median(), "p05": om.quantile(0.05), "p95": om.quantile(0.95), "sd": om.std()})
    st.index.name = "metric"
    st["n_members"] = obs.shape[1]
    st.to_csv(OBS_DIR / "metric_summary.csv")


def main():
    GD.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(exist_ok=True)
    if "--obs" in sys.argv:
        vendor_obs(sys.argv[sys.argv.index("--obs") + 1])
    pd.read_csv(OBS_DIR / "annual_summary.csv", index_col=0).to_csv(OUT / "obs_gsat_ensemble_summary.csv")
    st = pd.read_csv(OBS_DIR / "metric_summary.csv", index_col=0)
    print(f"obs 2011–2020 vs 1850–1900: {st.loc['dT_2011_2020', 'median']:.3f} [{st.loc['dT_2011_2020', 'p05']:.3f}, "
          f"{st.loc['dT_2011_2020', 'p95']:.3f}] (AR6: {AR6_2011_2020[0]} [{AR6_2011_2020[1]}, {AR6_2011_2020[2]}])")
    hist_rows = [{"model": "GSAT ensemble (obs)", "track": "obs", "n_members": int(st.n_members.iloc[0]),
                  **st["median"].to_dict(), **{f"{k}_p05": v for k, v in st.p05.items()},
                  **{f"{k}_p95": v for k, v in st.p95.items()}, **{f"{k}_sd": v for k, v in st.sd.items()}}]
    proj_rows, splice_all, checks_all, drift_all = [], [], [], []
    for model in MODELS:
        pi_cache = {}
        res, sp, ch, dr, ivs = build_model(model, pi_cache)
        splice_all += sp
        checks_all += ch
        drift_all += dr
        for exp, (df, used) in res.items():
            df.index.name = "year"
            df.to_csv(GD / f"{model}.{exp}.csv")
            if exp in SCEN:
                n = df.shape[1]
                # ensemble value = mean of per-member metrics over members that cover the
                # window (a skipna ensemble-mean series would change composition mid-window)
                mm = pd.DataFrame([metrics(df[c]) for c in df.columns])
                row = {"model": model, "track": exp, "n_members": n,
                       "splice": ",".join(sorted(set(used.values()))), **mm.mean().to_dict(),
                       **{f"{k}_n": int(mm[k].count()) for k in mm},
                       **{f"{k}_min": mm[k].min() for k in mm}, **{f"{k}_max": mm[k].max() for k in mm}}
                # Consistency with obs: model ensemble mean vs one obs realisation.
                # σ_iv from member spread if ≥5 members with the metric, else detrended piControl.
                for k in METRICS:
                    vals = mm[k].dropna()
                    if len(vals) >= 5:
                        s_iv, src = vals.std(ddof=1), "members"
                    else:
                        s_iv, src = iv_sigma(k, ivs.get(exp)), "piControl"
                    ne = max(len(vals), 1)
                    s = np.sqrt(st.sd[k] ** 2 + s_iv ** 2 * (1 + 1 / ne))
                    row[f"{k}_sigma_iv"], row[f"{k}_sigma_iv_src"] = s_iv, src
                    row[f"{k}_sigma"] = s
                    row[f"{k}_minus_obs"] = row[k] - st["median"][k]
                    row[f"{k}_z"] = row[f"{k}_minus_obs"] / s
                    # Fractional reduction of model warming that would make it consistent
                    # (z = 1.645): a yardstick for any GSAT-vs-GMST blending effect, which is
                    # not applied (project convention, METHODS §1.2).
                    if row[f"{k}_z"] > 1.645:
                        row[f"{k}_reduction_to_consistent"] = 1 - (st["median"][k] + 1.645 * s) / row[k]
                # Drift sensitivities (not the headline), z against the same σ:
                #   _driftcorr     remove the full-record control slope (the drift itself)
                #   _driftcorr_seg remove each member's contemporaneous-segment slope (full-record
                #                  where the segment is incomplete); mixes drift with centennial
                #                  variability, which the member spread already samples.
                dmap = {d["member"]: d for d in dr if d["track"] == exp and "control_trend_full_per_dec" in d}
                if dmap:
                    full = pd.Series({m: d["control_trend_full_per_dec"] for m, d in dmap.items()})
                    seg = pd.Series({m: (d["control_trend_segment_per_dec"]
                                         if np.isfinite(d["control_trend_segment_per_dec"])
                                         else d["control_trend_full_per_dec"]) for m, d in dmap.items()})
                    for k, (kind, y0, y1) in METRICS.items():
                        lever = (np.mean([y0, y1]) - np.mean(BASE)) / 10 if kind == "mean" else 1
                        for suf, sl in (("", full), ("_seg", seg)):
                            v = (mm[k].set_axis(df.columns) - (sl * lever).reindex(df.columns)).mean()
                            row[f"{k}_driftcorr{suf}"] = v
                            row[f"{k}_driftcorr{suf}_z"] = (v - st["median"][k]) / row[f"{k}_sigma"]
                hist_rows.append(row)
            else:
                r = {"model": model, "scenario": exp, "n_members": df.shape[1], "last_year": int(df.index.max())}
                for (a, b) in [(2041, 2060), (2081, 2100)]:
                    per = [metric(df[c], "mean", a, b) for c in df.columns]
                    r[f"dT_{a}_{b}"] = np.nanmean(per)
                    r[f"dT_{a}_{b}_min"], r[f"dT_{a}_{b}_max"] = np.nanmin(per), np.nanmax(per)
                proj_rows.append(r)
    pd.DataFrame(hist_rows).to_csv(OUT / "historical_metrics.csv", index=False)
    pd.DataFrame(proj_rows).to_csv(OUT / "projection_metrics.csv", index=False)
    pd.DataFrame(splice_all).to_csv(OUT / "splice_sensitivity.csv", index=False)
    ck = pd.DataFrame(checks_all)
    ck.to_csv(OUT / "splice_check.csv", index=False)
    dr = pd.DataFrame(drift_all)
    dr.to_csv(OUT / "drift_check.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 200):
        h = pd.DataFrame(hist_rows)
        print(h[["model", "track", "n_members", "splice", "dT_2012_2021", "dT_2012_2021_z", "dT_2015_2024",
                 "trend_1980_2021", "trend_1980_2021_z", "trend_1950_2021", "trend_1950_2021_z"]].round(3).to_string())
        print(pd.DataFrame(proj_rows).round(2).to_string())
        print("\nscenario runs not used:\n", ck[~ck.used].to_string())
        print("\nlargest |2021→2022 jump|:\n", ck.reindex(ck.jump_2021_2022.abs().sort_values().index[-5:]).to_string())
        print("\ndrift (°C added to 2012–2021 vs 1850–1900):\n",
              dr.groupby(["model", "track"]).agg(n=("member", "size"), drift_mean=("drift_in_dT_2012_2021", "mean"),
                                                 drift_min=("drift_in_dT_2012_2021", "min"),
                                                 drift_max=("drift_in_dT_2012_2021", "max"),
                                                 drift_full=("drift_in_dT_2012_2021_fullrecord", "first"),
                                                 seg=("segment", "first")).round(3).to_string())


if __name__ == "__main__":
    main()
