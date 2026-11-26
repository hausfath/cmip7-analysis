"""Earth's energy imbalance (EEI): CMIP7 models vs observational estimates (METHODS.md §2.9).

Model EEI = net TOA flux N = rsdt - rsut - rlut (global, annual, W m-2), with the parent control's
full-record OLS line removed at each member's contemporaneous control years (most of that is a
constant energy leak, not drift). historical -> piControl, esm-hist -> esm-piControl. 2022-2024 come
from the scenario joins used for GSAT (outputs/splice_check.csv, preference m > h > vl).
Models without published control radiation are not analysed.

Observations (data/obs/eei/, built by scripts/obs_eei.py), all as W m-2 of Earth's surface:
  wu        Wu et al. 2025: OLS rate of annual 0-2000 m OHC per member (72), + deep, / f_ocean
  miniere   Miniere et al. 2026: OLS rate of annual 0-2000 m OHC per product, + deep, / f_ocean
  hakuba    Hakuba et al. 2021: published EEI, hydrographic and geodetic methods (2005-2019 only)
  ceres     CERES EBAF Ed4.2.1: annual net TOA flux. Its level is tied to an in-situ ocean heating
            value for July 2005-June 2015 (Loeb et al. 2018), so it is shown but not used to test
            levels; its change between windows is independent of that anchor and is tested.

Estimators: a model is compared with OHC-based estimates through the OLS slope of its cumulative
heat (the same estimator applied to an OHC series) and with CERES through the window mean of N.

Outputs: derived/eei/<model>.<track>.csv (annual corrected N per member)
         outputs/eei_metrics.csv, eei_obs.csv, eei_tests.csv, eei_sensitivity.csv, eei_control.csv
"""
import json, pathlib, sys
import numpy as np
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # -I drops the script dir
import cmiplib as L  # noqa: E402

OUT = L.ROOT / "outputs"
DER = L.ROOT / "derived" / "eei"
OBS = L.ROOT / "data" / "obs" / "eei"
SEC = 365.25 * 86400
A_EARTH = 4 * np.pi * 6.371e6 ** 2
HIST_END, LAST = 2021, 2024
CONTROL = {"historical": "piControl", "esm-hist": "esm-piControl"}
SCEN = {"historical": ["scen7-m", "scen7-h", "scen7-vl"], "esm-hist": ["esm-scen7-m", "esm-scen7-h", "esm-scen7-vl"]}
# Ocean share of Earth's heat gain (von Schuckmann et al. 2023: "about 89 %" for 1971-2020 and
# 2006-2020); AR6 SPM A.4.2 "91%" as a sensitivity
F_OCEAN, F_OCEAN_ALT = 0.89, 0.91
# Ocean below 2000 m (von Schuckmann et al. 2023): 0.06 +/- 0.03 W m-2 from 1992, zero before.
# The +/- is treated as 1 sigma (the wider reading, which makes the test less likely to reject).
DEEP, DEEP_SD, DEEP_FROM = 0.06, 0.03, 1992
# Fixed, pre-registered windows. level: mean (CERES) or OLS rate of cumulative heat (OHC-based);
# change: mean 2015-2024 minus mean 2005-2014.
WINDOWS = {"2005_2019": (2005, 2019), "2015_2024": (2015, 2024)}
CHANGE = ((2005, 2014), (2015, 2024))
Z = 1.645


# ---------------------------------------------------------------- estimators
def ols_rate(heat, y0, y1):
    """OLS slope of an annual heat series (J m-2) over y0..y1 -> W m-2."""
    s = heat.loc[y0:y1].dropna()
    if len(s) != y1 - y0 + 1:
        return np.nan
    return np.polyfit(s.index.to_numpy(float), s.to_numpy(), 1)[0] / SEC


def cum_heat(n):
    """Annual-mean cumulative heat (J m-2) implied by annual-mean flux n: integral to mid-year."""
    return (n.cumsum() - n / 2) * SEC


def model_metrics(n):
    r = {}
    for k, (y0, y1) in WINDOWS.items():
        w = n.loc[y0:y1]
        full = w.notna().sum() == y1 - y0 + 1
        r[f"mean_{k}"] = w.mean() if full else np.nan
        r[f"rate_{k}"] = ols_rate(cum_heat(n.dropna()), y0, y1) if full else np.nan
    a, b = (n.loc[y0:y1] for y0, y1 in CHANGE)
    r["change"] = b.mean() - a.mean() if b.notna().sum() == 10 and a.notna().sum() == 10 else np.nan
    return r


def deep(y0, y1):
    return DEEP * np.mean([y >= DEEP_FROM for y in range(y0, y1 + 1)])


# ---------------------------------------------------------------- observations
def obs_estimates(mask="no_mask", f=F_OCEAN):
    rows = []
    wu = pd.read_csv(OBS / "wu2025_gohc_0-2000m_ZJ.csv", index_col=0) * 1e21 / A_EARTH
    mi = pd.read_csv(OBS / f"miniere2026_gohc_0-2000m_J_{mask}.csv", index_col=0).drop(columns="ensemble_mean") / A_EARTH
    for k, (y0, y1) in WINDOWS.items():
        for name, df in [("wu", wu), ("miniere", mi)]:
            r = pd.Series({c: ols_rate(df[c], y0, y1) for c in df}).dropna()
            if not len(r):
                continue
            e = (r + deep(y0, y1)) / f
            sd = np.sqrt(e.std(ddof=1) ** 2 + (DEEP_SD * deep(y0, y1) / DEEP / f) ** 2)
            rows.append({"source": name, "metric": f"rate_{k}", "value": e.mean(), "sd": sd,
                         "p05": e.quantile(0.05), "p95": e.quantile(0.95), "n": len(e),
                         "members": ",".join(r.index)})
    hk = pd.read_csv(OBS / "hakuba2021_eei_windows.csv")
    for r in hk[hk.window == "2005-2019"].itertuples():
        rows.append({"source": f"hakuba_{r.method}", "metric": "rate_2005_2019", "value": r.eei_mean,
                     "sd": (r.eei_p95 - r.eei_p05) / (2 * Z), "p05": r.eei_p05, "p95": r.eei_p95, "n": 1})
    ce = pd.read_csv(OBS / "ceres_ebaf_ed4.2.1_annual.csv", index_col=0).net_toa_Wm2
    for k, (y0, y1) in WINDOWS.items():
        rows.append({"source": "ceres", "metric": f"mean_{k}", "value": ce.loc[y0:y1].mean(), "sd": np.nan, "n": 1})
    (a0, a1), (b0, b1) = CHANGE
    rows.append({"source": "ceres", "metric": "change", "value": ce.loc[b0:b1].mean() - ce.loc[a0:a1].mean(),
                 "sd": 0.0, "n": 1})
    ch = ((mi.apply(lambda c: ols_rate(c, b0, b1)) - mi.apply(lambda c: ols_rate(c, a0, a1))) / f).dropna()
    rows.append({"source": "miniere", "metric": "change", "value": ch.mean(), "sd": ch.std(ddof=1),
                 "p05": ch.quantile(0.05), "p95": ch.quantile(0.95), "n": len(ch), "members": ",".join(ch.index)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- models
def net(model, exp, mem):
    sd, a = L.load_cmip7(model, exp, mem, "rsdt")
    su, _ = L.load_cmip7(model, exp, mem, "rsut")
    lu, _ = L.load_cmip7(model, exp, mem, "rlut")
    return (sd - su - lu).dropna(), a


def _complete():
    """Datasets whose every manifest file has been reduced (partially acquired ones are skipped)."""
    m = pd.read_csv(L.ROOT / "inventory" / "manifest_latest.csv")
    m["done"] = [(L.MONTHLY / k / (f[:-3] + ".csv")).exists() for k, f in zip(m.dskey, m.filename)]
    return set(m.groupby("dskey").done.all().pipe(lambda s: s[s].index))


COMPLETE = None


def has_rad(model, exp, mem):
    global COMPLETE
    if COMPLETE is None:
        COMPLETE = _complete()
    for v in ("rsdt", "rsut", "rlut"):
        try:
            if L.dskey(model, exp, mem, v) not in COMPLETE:
                return False
        except FileNotFoundError:
            return False
    return True


def control(model, cexp, want=None):
    mems = [m for m in L.members_cmip7(model, cexp, "rsdt") if has_rad(model, cexp, m)]
    if not mems:
        return None
    m = want if want in mems else mems[0]
    n, a = net(model, cexp, m)
    t = n.index.to_numpy(float)
    return {"member": m, "N": n, "attrs": a, "lin": np.polyfit(t, n.to_numpy(), 1),
            "quad": np.polyfit(t, n.to_numpy(), 2), "mean": float(n.mean())}


def ctl_variability(n, fit):
    """sigma of each metric from the control's residual about its fit (overlapping windows)."""
    x = n - np.polyval(fit, n.index.to_numpy(float))
    x.index = np.arange(len(x))
    out = {}
    for k, (y0, y1) in WINDOWS.items():
        L_ = y1 - y0 + 1
        out[f"mean_{k}"] = float(x.rolling(L_).mean().std())
        h = cum_heat(x)
        tt = np.arange(L_, dtype=float) - (L_ - 1) / 2
        w = tt / (tt ** 2).sum()
        out[f"rate_{k}"] = float(np.std(np.convolve(h.to_numpy(), w[::-1], mode="valid") / SEC, ddof=1))
    m10 = x.rolling(10).mean()
    out["change"] = float((m10 - m10.shift(10)).std())
    return out


def splice_parents():
    sc = pd.read_csv(OUT / "splice_check.csv")
    sc = sc[sc.used]
    return {(r.model, r.parent_member, r.scenario): r.member for r in sc.itertuples()}


def build(model, joins):
    series, ctl_rows, sens = {}, [], []
    for hexp, cexp in CONTROL.items():
        mems = [m for m in L.members_cmip7(model, hexp, "rsdt") if has_rad(model, hexp, m)]
        if not mems:
            continue
        first_attrs = net(model, hexp, mems[0])[1]
        c = control(model, cexp, first_attrs.get("parent_variant_label"))
        if c is None:
            ctl_rows.append({"model": model, "track": hexp, "control": cexp, "note": "no control radiation published"})
            continue
        cols, raw_cols = {}, {}
        for mem in mems:
            h, ha = net(model, hexp, mem)
            h = h.loc[:HIST_END]
            for s in SCEN[hexp]:
                smem = joins.get((model, mem, s))
                if smem and has_rad(model, s, smem):
                    x, _ = net(model, s, smem)
                    h = pd.concat([h, x.loc[HIST_END + 1:LAST]])
                    break
            by, off, _ = L.branch_year(ha, c["attrs"])
            cy = by + (h.index.to_numpy(float) - h.index.min())       # contemporaneous control years
            raw_cols[mem] = h
            cols[mem] = h - np.polyval(c["lin"], cy)
            for variant, ref in [("control_mean", np.full(len(cy), c["mean"])),
                                 ("control_quadratic", np.polyval(c["quad"], cy))]:
                sens.append({"model": model, "track": hexp, "member": mem, "variant": variant,
                             **model_metrics(h - ref)})
            ctl_rows.append({"model": model, "track": hexp, "member": mem, "control": cexp,
                             "control_member": c["member"], "control_years": f"{int(c['N'].index.min())}-{int(c['N'].index.max())}",
                             "control_mean_N": c["mean"], "control_slope_per_century": 100 * c["lin"][0],
                             "branch_year": by, "branch_offset_days": off,
                             "correction_2005_2019": float(np.polyval(c["lin"], cy[(h.index >= 2005) & (h.index <= 2019)]).mean())})
        df = pd.DataFrame(cols)
        series[hexp] = (df, ctl_variability(c["N"], c["lin"]), c)
    return series, ctl_rows, sens


def main():
    DER.mkdir(parents=True, exist_ok=True)
    joins = splice_parents()
    obs = obs_estimates()
    rows, tests, ctl_all, sens_all = [], [], [], []
    models = sorted({p.name.split(".")[4] for p in L.MONTHLY.glob("*.historical.*.rsdt.*")} |
                    {p.name.split(".")[4] for p in L.MONTHLY.glob("*.esm-hist.*.rsdt.*")})
    for model in models:
        series, ctl_rows, sens = build(model, joins)
        ctl_all += ctl_rows
        sens_all += sens
        for hexp, (df, iv, c) in series.items():
            df.round(5).to_csv(DER / f"{model}.{hexp}.csv", index_label="year")
            mm = pd.DataFrame([model_metrics(df[m]) for m in df])
            row = {"model": model, "track": hexp, "n_members": df.shape[1], "control": CONTROL[hexp],
                   "control_member": c["member"]}
            for k in mm:
                v = mm[k].dropna()
                row.update({f"{k}": v.mean() if len(v) else np.nan, f"{k}_min": v.min() if len(v) else np.nan,
                            f"{k}_max": v.max() if len(v) else np.nan, f"{k}_n": len(v), f"{k}_sigma_iv": iv[k]})
                if len(v) >= 5:
                    row[f"{k}_sigma_members"] = v.std(ddof=1)
            rows.append(row)
            for o in obs.itertuples():
                if o.metric not in mm or np.isnan(o.sd) or not row[f"{o.metric}_n"]:
                    continue
                n = row[f"{o.metric}_n"]
                sig = np.sqrt(o.sd ** 2 + iv[o.metric] ** 2 * (1 + 1 / n))
                tests.append({"model": model, "track": hexp, "metric": o.metric, "obs": o.source,
                              "model_value": row[o.metric], "obs_value": o.value, "sigma": sig,
                              "z": (row[o.metric] - o.value) / sig, "n": n})
    met = pd.DataFrame(rows)
    met.to_csv(OUT / "eei_metrics.csv", index=False)
    obs.to_csv(OUT / "eei_obs.csv", index=False)
    pd.DataFrame(ctl_all).to_csv(OUT / "eei_control.csv", index=False)
    t = pd.DataFrame(tests)
    # Verdict per model/track/metric: inconsistent only if |z| > 1.645 against EVERY observational
    # estimate available for that metric (same sign); otherwise consistent with at least one.
    def verdict(g):
        if (g.z > Z).all():
            return "higher than all tested"
        if (g.z < -Z).all():
            return "lower than all tested"
        return "consistent with at least one"
    v = t.groupby(["model", "track", "metric"]).apply(verdict, include_groups=False).rename("verdict").reset_index()
    t = t.merge(v, on=["model", "track", "metric"])
    t.to_csv(OUT / "eei_tests.csv", index=False)

    # sensitivities: control-correction variants, f_ocean, Miniere mask
    sens = pd.DataFrame(sens_all).groupby(["model", "track", "variant"]).mean(numeric_only=True).reset_index()
    alt = []
    for f_, mask, lab in [(F_OCEAN_ALT, "no_mask", "f_ocean_0.91"), (F_OCEAN, "GCOS_surface", "miniere_GCOS_mask")]:
        o = obs_estimates(mask, f_)
        alt.append(o.assign(variant=lab))
    pd.concat([sens.assign(kind="model"), pd.concat(alt).assign(kind="obs")]).to_csv(OUT / "eei_sensitivity.csv", index=False)

    pd.set_option("display.width", 200)
    print(obs[["source", "metric", "value", "sd", "p05", "p95", "n"]].round(3).to_string(index=False))
    print(met[["model", "track", "n_members", "rate_2005_2019", "mean_2005_2019", "rate_2015_2024",
               "mean_2015_2024", "change"]].round(3).to_string(index=False))
    print(v.to_string(index=False))


if __name__ == "__main__":
    main()
