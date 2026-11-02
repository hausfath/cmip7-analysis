"""Build the HTML report from outputs/*.csv and figures/*.png. Every number on the page is
read from the computed outputs; none is typed in.

  python -I scripts/build_report.py             ->  docs/index.html (+ docs/figures/*.png), the GitHub Pages site
  python -I scripts/build_report.py --artifact  ->  report/index.html, the retired claude.ai copy with a pointer
                                                    to the Pages site (report/ is not tracked)
"""
import html, json, pathlib, re, shutil, sys
import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
ARTIFACT = "--artifact" in sys.argv
OUT, FIG = ROOT / "outputs", ROOT / "figures"
REP = ROOT / ("report" if ARTIFACT else "docs")
REPO = "https://github.com/hausfath/cmip7-analysis"
SITE = "https://hausfath.github.io/cmip7-analysis/"
SNAPSHOT = "2026-10-24"
AR6_ECS_LIK, AR6_ECS_VL, AR6_TCR_LIK, AR6_TCR_VL = (2.5, 4.0), (2.0, 5.0), (1.4, 2.2), (1.2, 2.4)
TRACK = {"historical": "Concentration-driven", "esm-hist": "Emissions-driven"}


def f2(x, d=2):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{d}f}"


def esc(s):
    return html.escape(str(s))


def table(head, rows, cls=""):
    th = "".join(f"<th>{h}</th>" for h in head)
    tr = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'


def figure(name, caption, narrow=False):
    return (f'<figure{" class=narrow" if narrow else ""}><div class="plate"><img src="figures/{name}.png" alt="{esc(caption)}" loading="lazy"></div>'
            f"<figcaption>{caption}</figcaption></figure>")


def main():
    REP.mkdir(exist_ok=True)
    (REP / "figures").mkdir(exist_ok=True)

    man = pd.read_csv(ROOT / "inventory" / "manifest_latest.csv")
    dst = pd.read_csv(OUT / "cmip7_dataset_table.csv")
    c6 = pd.read_csv(OUT / "cmip6_sensitivity.csv")
    c7 = pd.read_csv(OUT / "cmip7_sensitivity.csv")
    lp = pd.read_csv(OUT / "lineage_pairs.csv")
    cs1 = lp.set_index("model").loc["CanESM5-1"]
    hm = pd.read_csv(OUT / "historical_metrics.csv")
    pm = pd.read_csv(OUT / "projection_metrics.csv")
    sp = pd.read_csv(OUT / "splice_sensitivity.csv")
    ref = pd.read_csv(ROOT / "reference" / "hausfather2022_cmip6_ecs_tcr.csv", comment="#")
    ref = ref[~ref.copy_of_other_model]

    n_models = man.source_id.nunique()
    n_files, gb = len(man), man["size"].sum() / 1e9
    n_stale = int(dst.files_ok_stale_index.sum())

    # ---- validation
    val = c6[c6.validation_set]
    dE = val.dECS_vs_zelinka.abs().max()
    dF = val.dF2x_vs_zelinka.abs().max()
    dL = val.dlambda_vs_zelinka.abs().max()
    dT_ar6 = val.dTCR_linear_vs_AR6.abs().max()
    dT_h22 = val.dTCR_linear_vs_H22.abs().max()
    dT_rm = val.dTCR_rm21_vs_H22.abs().max()
    allE = c6.dECS_vs_zelinka.abs().max()

    # ---- sensitivity
    lp = lp.sort_values("ECS", ascending=False)
    n7 = len(lp)
    n_lower = int((lp.ECS_change == "lower").sum())
    n_lower1 = int((lp.ECS_change_single_range == "lower").sum())
    n_untested = int((lp.ECS_change == "untested").sum())
    n_above_lik = int((lp.ECS > AR6_ECS_LIK[1]).sum())
    n_above_vl = int((lp.ECS > AR6_ECS_VL[1]).sum())
    tcr = lp.dropna(subset=["TCR"])
    n_tcr_above_lik = int((tcr.TCR > AR6_TCR_LIK[1]).sum())
    e6 = ref.ECS150.dropna()
    t6 = ref.TCR.dropna()
    frac6_above_vl = (e6 > AR6_ECS_VL[1]).mean()
    frac6_above_lik = (e6 > AR6_ECS_LIK[1]).mean()

    # ---- historical
    obs = hm[hm.track == "obs"].iloc[0]
    mh = hm[hm.track != "obs"]
    t = mh.dropna(subset=["dT_2012_2021_z"]).groupby("model").dT_2012_2021_z
    n_tested = t.ngroups
    n_hot = int((t.min() > 1.645).sum())
    hot_models = sorted(t.min()[t.min() > 1.645].index)
    split = sorted(t.min()[(t.max() > 1.645) & (t.min() <= 1.645)].index)
    cons = sorted(t.max()[t.max().abs() <= 1.645].index)
    untested = sorted(set(mh.model) - set(t.groups))
    ts = mh.dropna(subset=["dT_2012_2021_driftcorr_seg_z"]).groupby("model").dT_2012_2021_driftcorr_seg_z
    n_hot_seg = int((ts.min() > 1.645).sum())
    drop_seg = sorted(set(hot_models) - set(ts.min()[ts.min() > 1.645].index))
    tf = mh.dropna(subset=["dT_2012_2021_driftcorr_z"]).groupby("model").dT_2012_2021_driftcorr_z
    n_hot_full = int((tf.min() > 1.645).sum())
    tt = mh.dropna(subset=["trend_1980_2021_z"]).groupby("model").trend_1980_2021_z
    n_trend_hot = int((tt.min() > 1.645).sum())

    # ---- projections
    proj = pm.copy()
    p_h = proj[proj.scenario.str.endswith("scen7-h")]
    p_vl = proj[proj.scenario.str.endswith("scen7-vl")]
    n_proj = proj.model.nunique()
    warm_any = set(t.max()[t.max() > 1.645].index)
    n_proj_warm = len(set(proj.model) & warm_any)
    n_late_larger = int((lp.ECS_21_150 > lp.ECS).sum())
    inst = man.drop_duplicates("source_id").institution_id
    n_ukcaau = int(inst.isin(["CCCma", "UKNCSP", "MOHC", "NERC", "ACCESS-Consortium", "CSIRO"]).sum())
    sp_rng = sp.groupby(["model", "track"]).dT_2015_2024.agg(lambda x: x.max() - x.min()).max()

    # ---------------------------------------------------------------- tables
    sens_rows = []
    for _, r in lp.iterrows():
        c = c7[c7.model == r.model].iloc[0]
        sens_rows.append([
            f"<b>{esc(r.model)}</b>", esc(r.member),
            f"{f2(r.ECS)} <span class='rng'>{f2(r.ECS_p05)}–{f2(r.ECS_p95)}</span>",
            f2(r.ECS_21_150), f2(r.F2x), f2(r["lambda"]), f2(r.TCR),
            f"{esc(r.predecessor)}", f2(r.pred_ECS), f2(r.pred_TCR),
            f"<span class='chip {'lo' if r.ECS_change == 'lower' else ('hi' if r.ECS_change == 'higher' else 'nc')}'>"
            f"{f2(r.dECS, 2) if np.isfinite(r.dECS) else '–'}</span>"])
    sens_tab = table(["CMIP7 model", "Member", "ECS (5–95%)", "ECS yr 21–150", "ERF₂ₓ", "λ", "TCR",
                      "CMIP6 predecessor", "ECS", "TCR", "ΔECS"], sens_rows, "num")

    val_rows = []
    for _, r in c6.iterrows():
        val_rows.append([f"<b>{esc(r.model)}</b> {esc(r.member)}" + (" <span class='tag'>validation</span>"
                         if r.validation_set else ""),
                         f2(r.ECS_rm21, 3), f2(r.ECS_zelinka, 3), f2(r.dECS_vs_zelinka, 3),
                         f2(r.F2x_rm21, 3), f2(r.dF2x_vs_zelinka, 3),
                         f2(r.TCR_linear, 3), f2(r.TCR_AR6, 2), f2(r.dTCR_linear_vs_AR6, 3),
                         f2(r.TCR_rm21, 3), f2(r.dTCR_rm21_vs_H22, 3)])
    val_tab = table(["CMIP6 run", "ECS (ours)", "ECS Zelinka v2.2", "Δ", "ERF₂ₓ", "Δ", "TCR linear",
                     "TCR AR6", "Δ", "TCR rm21", "Δ vs H22"], val_rows, "num")

    hist_rows = [["<b>GMST ensemble</b>", "obs", f"{int(obs.n_members):,}",
                  f"{f2(obs.dT_2012_2021)} <span class='rng'>{f2(obs.dT_2012_2021_p05)}–{f2(obs.dT_2012_2021_p95)}</span>",
                  "", f"{f2(obs.dT_2015_2024)}", f"{f2(obs.trend_1980_2021, 3)}", "", f"{f2(obs.trend_1950_2021, 3)}",
                  ""]]
    for _, r in mh.sort_values("dT_2012_2021", ascending=False).iterrows():
        z = r.dT_2012_2021_z
        zc = "hi" if z > 1.645 else ("lo" if z < -1.645 else "nc")
        zs = f"<span class='chip {zc}'>{f2(z, 1)}</span>" if np.isfinite(z) else "<span class='muted'>n/a</span>"
        zt = r.trend_1980_2021_z
        ztc = "hi" if zt > 1.645 else ("lo" if zt < -1.645 else "nc")
        zts = f"<span class='chip {ztc}'>{f2(zt, 1)}</span>" if np.isfinite(zt) else "<span class='muted'>n/a</span>"
        hist_rows.append([f"<b>{esc(r.model)}</b>", TRACK[r.track], int(r.n_members),
                          f"{f2(r.dT_2012_2021)} <span class='rng'>{f2(r.dT_2012_2021_min)}–{f2(r.dT_2012_2021_max)}</span>",
                          zs, f2(r.dT_2015_2024), f2(r.trend_1980_2021, 3), zts, f2(r.trend_1950_2021, 3),
                          f2(r.get("dT_2012_2021_driftcorr", np.nan))])
    hist_tab = table(["Model", "Track", "n", "2012–2021 (°C, range)", "z", "2015–2024 (°C)",
                      "Trend 1980–2021 (°C/dec)", "z", "Trend 1950–2021", "2012–2021 drift-corr."], hist_rows, "num")

    proj_rows = []
    for _, r in proj.sort_values(["scenario", "dT_2081_2100"], ascending=[True, False]).iterrows():
        proj_rows.append([f"<b>{esc(r.model)}</b>", esc(r.scenario), int(r.n_members),
                          f"{f2(r.dT_2041_2060)} <span class='rng'>{f2(r.dT_2041_2060_min)}–{f2(r.dT_2041_2060_max)}</span>",
                          f"{f2(r.dT_2081_2100)} <span class='rng'>{f2(r.dT_2081_2100_min)}–{f2(r.dT_2081_2100_max)}</span>"])
    proj_tab = table(["Model", "Scenario", "n", "2041–2060 (°C)", "2081–2100 (°C)"], proj_rows, "num")

    hold = []
    for src, g in dst.groupby("source_id"):
        def n(e, var="tas"):
            x = g[(g.experiment_id == e) & g.variable.str.startswith(var + "_")]
            return x.variant_label.nunique()
        sc = "/".join(str(n(f"scen7-{s}") or "–") for s in ("h", "m", "vl"))
        esc_ = "/".join(str(n(f"esm-scen7-{s}") or "–") for s in ("h", "m", "vl"))
        hold.append([f"<b>{esc(src)}</b>", esc(g.institution_id.iloc[0]), n("historical") or "–", n("esm-hist") or "–",
                     sc, esc_, "✓" if n("abrupt-4xCO2", "rlut") else "–", "✓" if n("1pctCO2") else "–"])
    hold_tab = table(["Model", "Institution", "historical", "esm-hist", "scen7 h/m/vl", "esm-scen7 h/m/vl",
                      "abrupt-4xCO2 radiation", "1pctCO2"], hold, "num")

    hot_txt = ", ".join(hot_models)
    head = ('<meta name="claude:viewport-fit-sides" content="cover">' if ARTIFACT else
            '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            '<meta name="description" content="CMIP7 historical warming and climate sensitivity compared with '
            'observations and CMIP6.">')
    moved = (f'<div class="moved">This page has moved to <a href="{SITE}">{SITE}</a>, with code and data at '
             f'<a href="{REPO}">{REPO}</a>. This copy is no longer updated.</div>' if ARTIFACT else "")
    page = f"""{head}
<title>CMIP7 First Look</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans+Condensed:wght@500;600&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: single reading column (~72ch) for prose; figures and tables break out to a wider measure. */
:root {{
  --bg: #f6f7f8; --surface: #ffffff; --ink: #14171a; --ink-2: #4a525a; --ink-3: #737b83;
  --rule: #dde1e5; --accent: #1f63b8; --hot: #b9471c; --hot-bg: #fbe9e1; --cool: #1d6b8f; --cool-bg: #e3f0f6;
  --nc-bg: #eceef0; --plate: #ffffff;
  --display: "IBM Plex Sans Condensed", "Arial Narrow", system-ui, sans-serif;
  --body: "IBM Plex Sans", system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, "SF Mono", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{ :root:not([data-theme="light"]) {{
  --bg: #111417; --surface: #181c20; --ink: #e8ebee; --ink-2: #b3bac1; --ink-3: #8a9299;
  --rule: #2b3137; --accent: #6fa8f0; --hot: #f08a5d; --hot-bg: #3a2219; --cool: #6cc0e0; --cool-bg: #16303b;
  --nc-bg: #252a2f; --plate: #f4f5f6; color-scheme: dark }} }}
:root[data-theme="dark"] {{
  --bg: #111417; --surface: #181c20; --ink: #e8ebee; --ink-2: #b3bac1; --ink-3: #8a9299;
  --rule: #2b3137; --accent: #6fa8f0; --hot: #f08a5d; --hot-bg: #3a2219; --cool: #6cc0e0; --cool-bg: #16303b;
  --nc-bg: #252a2f; --plate: #f4f5f6; color-scheme: dark }}
body {{ background: var(--bg); color: var(--ink); font: 16px/1.6 var(--body);
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px)); }}
main {{ max-width: 1080px; margin: 0 auto; padding-block: 40px 72px; display: grid;
  grid-template-columns: minmax(0, 1fr); gap: 56px; }}
section {{ display: grid; gap: 18px; min-width: 0; }}
.prose {{ max-width: 72ch; display: grid; gap: 12px; }}
.prose p, .prose li {{ margin: 0; }}
h1, h2, h3 {{ font-family: var(--display); font-weight: 600; text-wrap: balance; margin: 0; line-height: 1.15; }}
h1 {{ font-size: clamp(2rem, 4.5vw, 3rem); letter-spacing: -0.01em; }}
h2 {{ font-size: 1.6rem; padding-top: 18px; border-top: 1px solid var(--rule); }}
h3 {{ font-size: 1.15rem; }}
.eyebrow {{ font: 500 0.78rem/1 var(--mono); letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink-3); }}
.lede {{ font-size: 1.15rem; color: var(--ink-2); max-width: 66ch; margin: 0; }}
header {{ display: grid; gap: 14px; }}
.meta {{ display: flex; flex-wrap: wrap; gap: 8px 22px; font: 0.85rem/1.4 var(--mono); color: var(--ink-3); }}
.findings {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 420px), 1fr)); gap: 14px; }}
.finding {{ background: var(--surface); border: 1px solid var(--rule); border-radius: 6px; padding: 16px 18px;
  display: grid; gap: 6px; align-content: start; min-width: 0; }}
.finding .k {{ font: 600 1.7rem/1.1 var(--display); font-variant-numeric: tabular-nums; }}
.finding p {{ margin: 0; color: var(--ink-2); font-size: 0.95rem; }}
figure {{ margin: 0; display: grid; gap: 10px; min-width: 0; }}
figure.narrow {{ max-width: 680px; }}
.plate {{ background: var(--plate); border-radius: 6px; padding: 10px; border: 1px solid var(--rule); }}
.plate img {{ display: block; width: 100%; height: auto; }}
figcaption {{ font-size: 0.9rem; color: var(--ink-2); max-width: 90ch; }}
.tw {{ overflow-x: auto; border: 1px solid var(--rule); border-radius: 6px; background: var(--surface); }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.86rem; }}
th, td {{ padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--rule); white-space: nowrap; vertical-align: top; }}
th {{ font: 500 0.74rem/1.3 var(--mono); text-transform: uppercase; letter-spacing: 0.04em; color: var(--ink-3);
  background: var(--bg); white-space: normal; min-width: 60px; }}
tbody tr:last-child td {{ border-bottom: 0; }}
table.num td {{ font-variant-numeric: tabular-nums; }}
.rng {{ color: var(--ink-3); font-size: 0.8rem; }}
.muted {{ color: var(--ink-3); }}
.chip {{ display: inline-block; padding: 1px 7px; border-radius: 10px; font-weight: 500; background: var(--nc-bg); }}
.chip.hi {{ background: var(--hot-bg); color: var(--hot); }}
.chip.lo {{ background: var(--cool-bg); color: var(--cool); }}
.tag {{ font: 500 0.7rem/1 var(--mono); padding: 2px 6px; border-radius: 3px; background: var(--nc-bg); color: var(--ink-2);
  text-transform: uppercase; letter-spacing: 0.05em; }}
code {{ font: 0.88em var(--mono); background: var(--nc-bg); padding: 1px 4px; border-radius: 3px; }}
a {{ color: var(--accent); }}
a:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
ul {{ margin: 0; padding-left: 1.2em; display: grid; gap: 6px; }}
.note {{ font-size: 0.88rem; color: var(--ink-2); }}
.moved {{ border: 1px solid var(--accent); border-radius: 4px; padding: 10px 14px; margin-bottom: 20px; font-size: 0.95rem; }}
</style>
{"" if ARTIFACT else "</head>\n<body>"}
<main>
<header>
  {moved}
  <div class="eyebrow">Working note · ESGF snapshot {SNAPSHOT}</div>
  <h1>CMIP7 first look: warming and climate sensitivity</h1>
  <p class="lede">The first {n_models} CMIP7 models on ESGF, compared with the GMST ensemble over the historical
  period and with CMIP6 for effective climate sensitivity (ECS) and transient climate response (TCR). The sample is
  small ({n_ukcaau} of {n_models} models from UK, Canadian and Australian centres), so these are model-by-model results, not an ensemble
  assessment.</p>
  <div class="meta"><span>{n_models} models</span><span>{n_files:,} files · {gb:.0f} GB streamed</span>
  <span>{len(dst)} datasets · all CC-BY-4.0</span><span>temperatures in °C vs 1850–1900</span>
  <span><a href="{REPO}">Code and data</a></span></div>
</header>

<section>
  <h2>Key results</h2>
  <div class="findings">
    <div class="finding"><span class="eyebrow">ECS vs CMIP6 predecessor</span><span class="k">{n_lower} of {n7} lower</span>
      <p>{n_lower} of {n7} CMIP7 models have a lower 150-yr Gregory ECS than their nearest CMIP6 predecessor, with
      non-overlapping 5–95% regression ranges ({n_lower1} under the looser test that only the predecessor's value lies
      outside the CMIP7 range{f"; {n_untested} untested" if n_untested else ""}). CMIP7 ECS spans
      {f2(lp.ECS.min())}–{f2(lp.ECS.max())} °C.</p></div>
    <div class="finding"><span class="eyebrow">ECS vs AR6</span><span class="k">{n_above_lik} of {n7} above likely</span>
      <p>{n_above_lik} of {n7} remain above the AR6 likely range (2.5–4 °C) and {n_above_vl} above very likely (2–5 °C).
      For comparison, {frac6_above_lik:.0%} of CMIP6 models were above 4 °C and {frac6_above_vl:.0%} above 5 °C.</p></div>
    <div class="finding"><span class="eyebrow">Historical warming 2012–2021</span><span class="k">{n_hot} of {n_tested} too warm</span>
      <p>The GMST ensemble gives {f2(obs.dT_2012_2021)} °C. {n_hot} of {n_tested} testable models warm significantly
      more in every published track ({esc(hot_txt)}).{(' ' + esc(', '.join(split)) + ' does so in one track only.') if split else ''}
      Removing the control-run drift leaves {n_hot_full} of {tf.ngroups}; removing each member's contemporaneous
      control-run trend, which also contains centennial variability, gives {n_hot_seg} of {ts.ngroups}{(' (' + esc(', '.join(drop_seg)) + ' drops out)') if drop_seg else ''}.</p></div>
    <div class="finding"><span class="eyebrow">Code validation</span><span class="k">≤ {dE:.3f} °C</span>
      <p>The same code reproduces published CMIP6 ECS (Zelinka et al. 2020) to ≤ {dE:.3f} °C and TCR (AR6 Table 7.SM.5)
      to ≤ {dT_ar6:.3f} °C on the three validation models.</p></div>
  </div>
</section>

<section>
  <h2>Climate sensitivity</h2>
  <div class="prose">
    <p>ECS is the 150-year Gregory estimate from <code>abrupt-4xCO2</code>, with anomalies taken against a 21-yr running
    mean of the contemporaneous <code>piControl</code> (Zelinka et al. 2020 convention). TCR is the mean of years 61–80 of
    <code>1pctCO2</code>, with the control drift removed by a linear fit (the ESMValTool convention behind AR6 Table 7.SM.5).
    Each CMIP7 model is paired with its nearest CMIP6 predecessor by model family. The pairings are our judgement, not
    an official lineage, and the predecessors are recomputed with the same code so the differences are method-consistent.
    CanESM5-1 is CanESM5.1-p2 (<a href="https://doi.org/10.5194/gmd-16-6553-2023">Sigmond et al. 2023</a>), whose atmosphere
    was retuned to an effective climate sensitivity of 4.09 K (their value; {cs1.ECS:.2f} °C here) and which is
    also a CMIP6 model (CMIP6 source_id <code>CanESM5-1</code>, release year 2022). Its lower ECS against CanESM5 ({cs1.pred_ECS:.2f} °C) therefore dates from that retuning, not from CMIP7.</p>
    <p>Across CMIP6, {len(e6)} models had a median ECS of {e6.median():.2f} °C and {len(t6)} a median TCR of {t6.median():.2f} °C
    (Hausfather et al. 2022 compilation). {n_tcr_above_lik} of {len(tcr)} CMIP7 models with a 1pctCO2 run have TCR above
    the AR6 likely range (1.4–2.2 °C).</p>
  </div>
  {figure("fig04_ecs_tcr_vs_cmip6", "ECS and TCR for each CMIP7 model (blue) with its CMIP6 predecessor (hollow) and the CMIP6 distribution (grey). Shading: AR6 likely and very likely ranges; dashed line: AR6 best estimate.")}
  {sens_tab}
  <p class="note">ECS range: fixed-design circular block bootstrap of the regression residuals (10-yr blocks, 2,000
  resamples); it reflects regression sampling only. "ECS yr 21–150" is the fit excluding the first 20 years; it is larger
  for {n_late_larger} of {n7} models, the usual pattern-effect curvature. ΔECS is coloured only when the two models' 5–95% ranges do not
  overlap. The ranges come from one realisation each and exclude realisation-to-realisation spread in ECS, so they
  understate the true uncertainty of each difference.</p>
  {figure("fig05_gregory", "Gregory plots: annual ΔN against ΔT for years 1–150, CMIP7 (blue) and CMIP6 predecessor (grey).")}
  {figure("fig06_ecs_tcr_scatter", "ECS against TCR for CMIP7 models on the CMIP6 distribution.", narrow=True)}
</section>

<section>
  <h2>Historical warming</h2>
  <div class="prose">
    <p>Each member is expressed relative to its own 1850–1900 mean, and each model is summarised by the mean of its
    members' metrics. Concentration-driven (<code>historical</code>) and emissions-driven (<code>esm-hist</code>) runs are
    kept separate. Runs end in 2021; 2022–2024 come from each member's own scenario continuation, joined by the parent
    metadata in the scenario files. The choice of scenario changes ensemble-mean 2015–2024 warming by at most
    {sp_rng:.3f} °C.</p>
    <p><b>z</b> compares the model ensemble mean with the observations, given observational uncertainty (spread of the
    10,000-member GMST ensemble) and internal variability (member spread with ≥ 5 members, otherwise the control run).
    A model counts as warmer than observed if z &gt; 1.645 (one-sided 5%) in every published track. Models consistent with
    observations in all tracks: {esc(', '.join(cons))}. Not tested (one member, no control run): {esc(', '.join(untested))}.
    For the 1980–2021 trend, {n_trend_hot} of {tt.ngroups} testable models warm faster than observed
    ({f2(obs.trend_1980_2021, 3)} °C per decade).</p>
  </div>
  {figure("fig01_historical_small_multiples", "Annual tas for each CMIP7 model against the GMST ensemble (black, 5–95% grey).")}
  {figure("fig02_historical_metrics", "2012–2021 warming and 1980–2021 trend, model ensemble mean and member range.")}
  {hist_tab}
  <p class="note">No blending or masking adjustment is applied to model tas (project convention: the GMST ensemble is treated
  as equivalent to GSAT). Any such adjustment would lower model warming, so it can only shrink the warm biases shown here.
  <code>historical_metrics.csv</code> lists the fractional reduction in model warming that would make each warm model
  consistent. The drift-corrected column removes the full-record control-run trend (the drift itself). Removing each member's
  contemporaneous 172-yr control trend instead, which mixes drift with centennial variability, is a second sensitivity
  (<code>_driftcorr_seg</code> columns).</p>
</section>

<section>
  <h2>Projections</h2>
  <div class="prose"><p>{n_proj} models have scenario runs joined to a published historical parent so far, through 2100 for scen7-h,
  scen7-m and scen7-vl.
  Under scen7-h, 2081–2100 warming ranges from {f2(p_h.dT_2081_2100.min())} to {f2(p_h.dT_2081_2100.max())} °C across
  models and tracks; under scen7-vl, from {f2(p_vl.dT_2081_2100.min())} to {f2(p_vl.dT_2081_2100.max())} °C. {n_proj_warm} of
  these {n_proj} warm more than observed over 2012–2021 in at least one track, so they should not be read as a CMIP7
  projection range.</p></div>
  {figure("fig03_projections", "Ensemble-mean tas by scenario; solid lines concentration-driven, dashed emissions-driven.")}
  {proj_tab}
</section>

<section>
  <h2>Validation and data checks</h2>
  <div class="prose">
    <p>Before any CMIP7 numbers were computed, the code was run on three CMIP6 models with published values. The
    pre-registered tolerance was ±0.05 °C. ECS matched Zelinka et al. (2020) v2.2 to ≤ {dE:.3f} °C (ERF₂ₓ ≤ {dF:.3f} W m⁻²,
    λ ≤ {dL:.3f} W m⁻² K⁻¹); across all {len(c6)} recomputed CMIP6 runs the largest ECS difference is {allE:.3f} °C. TCR with
    the linear control fit matched AR6 Table 7.SM.5 to ≤ {dT_ar6:.3f} °C and Hausfather et al. (2022) to ≤ {dT_h22:.3f} °C. With a
    running-mean control instead, TCR was up to {dT_rm:.3f} °C high, so the linear convention is used for TCR.</p>
  </div>
  {val_tab}
  <div class="prose"><ul>
    <li>{n_files - n_stale:,} of {n_files:,} files matched their ESGF sha256 checksum. The other {n_stale} (CanESM5-1 piControl)
    were rewritten on the CCCma node after indexing; they are complete and readable but cannot be checksum-verified, so
    CanESM5-1's ECS rests on the server copy.</li>
    <li>Scenario runs are joined to their historical parent by metadata, not variant label: some UKESM1-3-LL esm-scen7 f1
    members branch from f2 esm-hist members. Runs whose declared parent is unpublished are not used
    (<code>splice_check.csv</code>).</li>
    <li>Version tags: CanESM5-1 datasets carry <code>v20190429</code>, but their global attributes identify CMIP7 output
    created in 2026. UKESM1-3-LL <code>historical</code> is forcing variant f2, while its control and idealised runs are f1.</li>
    <li>Some datasets in the ESGF inventory are daily-maximum or daily-minimum brandings (<code>tmaxavg</code>,
    <code>tmax</code>) rather than means, and are excluded. This removes an apparent CNRM-ESM2-1e <code>esm-hist</code> and
    CAS-ESM2-1-CSM <code>scen7-vl</code>.</li>
  </ul></div>
  {hold_tab}
</section>

<section>
  <h2>Caveats</h2>
  <div class="prose"><ul>
    <li>The sample is early and non-random. Counts and ranges here describe these {n_models} models, not CMIP7.</li>
    <li>150-yr Gregory ECS is an effective sensitivity and understates the long-run equilibrium for models with strong
    pattern effects.</li>
    <li>CMIP7 forcing differs from CMIP6, so CMIP6→CMIP7 differences in historical warming mix model and forcing changes.
    ECS and TCR come from idealised CO₂ runs and are not affected.</li>
    <li>Emissions-driven runs simulate their own CO₂, so differences between their warming and that of the
    concentration-driven runs partly reflect CO₂.</li>
    <li>The consistency test assumes realistic model internal variability. For models with 2–3 members it relies on the
    control run.</li>
  </ul></div>
</section>

<section>
  <h2>Sources</h2>
  <div class="prose note"><ul>
    <li>CMIP7 output via ESGF, per the <a href="https://wcrp-cmip.github.io/cmip7-guidance/docs/CMIP7/Guidance_for_users/">CMIP7 user guidance</a>.
    Dataset IDs, versions and tracking IDs: <code>outputs/cmip7_dataset_table.csv</code>.</li>
    <li>Hausfather, Marvel, Schmidt, Nielsen-Gammon &amp; Zelinka (2022), <a href="https://doi.org/10.1038/d41586-022-01192-2">Nature 605, 26–29</a>;
    data <a href="https://zenodo.org/records/6476375">Zenodo 6476375</a>.</li>
    <li>Zelinka et al. (2020), GRL; values v2.2 from <a href="https://github.com/mzelinka/cmip56_forcing_feedback_ecs">github.com/mzelinka/cmip56_forcing_feedback_ecs</a>.</li>
    <li>Sigmond et al. (2023), <a href="https://doi.org/10.5194/gmd-16-6553-2023">GMD 16, 6553</a>.</li>
    <li>IPCC AR6 WGI Chapter 7 and Table 7.SM.5; SPM A.1.2 for 2011–2020 warming.</li>
    <li>CMIP6 model output from the Pangeo <code>gs://cmip6</code> archive.</li>
  </ul>
  <p>Full methods, parameter choices and sensitivity tests: <a href="{REPO}/blob/main/METHODS.md">METHODS.md</a>.
  Code, derived data and the pipeline to reproduce every number and figure: <a href="{REPO}">{REPO}</a>.</p></div>
</section>
</main>
{"" if ARTIFACT else "</body>\n</html>"}
"""
    for name in sorted(set(re.findall(r'src="figures/([\w-]+\.png)"', page))):
        shutil.copy(FIG / name, REP / "figures" / name)
    (REP / "index.html").write_text(page)
    print("wrote", REP / "index.html")


if __name__ == "__main__":
    main()
