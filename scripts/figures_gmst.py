"""Warming figures from derived/gmst and outputs/*.csv (no raw data needed).

  python -I scripts/figures_gmst.py

fig01_historical_small_multiples  per-model GSAT vs the GSAT ensemble, 1850–2024
fig02_historical_metrics          2012–2021 warming and 1980–2021 trend per model × track
fig02b_historical_metrics_cmip6   the same metrics for CMIP6 (Hausfather et al. 2022 data, one member per
                                  model); descriptive only, no consistency test
fig03_projections                 ensemble-mean GSAT by scenario, 1850–2100
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import figstyle as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
GD, OUT, FIG = ROOT / "derived" / "gmst", ROOT / "outputs", ROOT / "figures"
MODELS = ["CanESM6-0-MR", "CanESM5-1", "UKESM1-3-LL", "UKCM2-0-LL", "UKCM2a-0-HH", "ACCESS-ESM1-6",
          "CAS-ESM2-1-CSM", "EC-Earth3-ESM-1-1", "CNRM-ESM2-1e"]
OBS_END = 2025


def load(model, exp):
    f = GD / f"{model}.{exp}.csv"
    return pd.read_csv(f, index_col=0) if f.exists() else None


def fig01(obs, hm):
    fig, axs = plt.subplots(3, 3, figsize=(S.WIDTH, 150 * S.MM), sharex=True, sharey=True)
    for ax, m in zip(axs.flat, MODELS):
        ax.fill_between(obs.index, obs.p05, obs.p95, color=S.OBS_BAND, lw=0, alpha=0.7)
        ax.plot(obs.index, obs["median"], color=S.OBS, lw=0.9)
        lines = []
        for exp in ("historical", "esm-hist"):
            d = load(m, exp)
            if d is None:
                continue
            d = d.loc[:2024]
            c = S.TRACK_COLOR[exp]
            if d.shape[1] > 1:
                ax.fill_between(d.index, d.min(axis=1), d.max(axis=1), color=c, alpha=0.18, lw=0)
            ax.plot(d.index, d.mean(axis=1), color=c, lw=1.0)
            r = hm[(hm.model == m) & (hm.track == exp)].iloc[0]
            lines.append(f"{ {'historical': 'Conc.', 'esm-hist': 'Emis.'}[exp]} {r.dT_2012_2021:.2f} °C (n={int(r.n_members)})")
        ax.set_title(m, loc="left", fontweight="bold", pad=3)
        ax.text(0.03, 0.95, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=6, color=S.INK2,
                linespacing=1.4)
        ax.set_xlim(1850, 2025)
        ax.set_ylim(-0.6, 2.2)
    for ax in axs[:, 0]:
        ax.set_ylabel("°C vs 1850–1900")
    o = hm[hm.track == "obs"].iloc[0]
    handles = [Line2D([], [], color=S.OBS, lw=1.2, label=f"GSAT ensemble median (2012–2021: {o.dT_2012_2021:.2f} °C)"),
               Patch(color=S.OBS_BAND, label="GSAT ensemble 5–95%"),
               Line2D([], [], color=S.TRACK_COLOR["historical"], lw=1.4, label="Concentration-driven (historical) mean"),
               Line2D([], [], color=S.TRACK_COLOR["esm-hist"], lw=1.4, label="Emissions-driven (esm-hist) mean"),
               Patch(color="#999", alpha=0.35, label="Member range")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.045))
    fig.tight_layout(rect=(0, 0.02, 1, 0.9))
    S.titles(fig, "CMIP7 models published so far vs observed global temperature",
             "Annual global surface air temperature (tas), °C relative to 1850–1900. Panel text: 2012–2021 mean "
             "warming.\n2022–2024 from each member's own scenario continuation (scen7-m preferred, else h, else vl) where published.", gap=0.04)
    S.save(fig, FIG / "fig01_historical_small_multiples")


def fig02(hm):
    obs = hm[hm.track == "obs"].iloc[0]
    rows = hm[hm.track != "obs"].copy()
    rows["label"] = rows.model + "  (" + rows.track.map({"historical": "conc.", "esm-hist": "emis."}) + \
        ", n=" + rows.n_members.astype(int).astype(str) + ")"
    rows = rows.sort_values("dT_2012_2021")
    y = np.arange(len(rows))
    fig, axs = plt.subplots(1, 2, figsize=(S.WIDTH, 95 * S.MM), sharey=True)
    for ax, k, unit, lab in [(axs[0], "dT_2012_2021", "°C", "2012–2021 warming (°C vs 1850–1900)"),
                             (axs[1], "trend_1980_2021", "°C/decade", "1980–2021 trend (°C per decade)")]:
        ax.axvspan(obs[f"{k}_p05"], obs[f"{k}_p95"], color=S.OBS_BAND, alpha=0.6, lw=0)
        ax.axvline(obs[k], color=S.OBS, lw=0.9)
        for yi, (_, r) in zip(y, rows.iterrows()):
            c = S.TRACK_COLOR[r.track]
            if r.n_members > 1:
                ax.plot([r[f"{k}_min"], r[f"{k}_max"]], [yi, yi], color=c, lw=1.6, alpha=0.45,
                        solid_capstyle="round")
            ax.plot(r[k], yi, "o", ms=4.5, color=c, mec="white", mew=0.6)
        ax.set_xlabel(lab)
        ax.grid(axis="y", visible=False)
        ax.text(obs[k], len(rows) - 0.3, f" GSAT ensemble {obs[k]:.2f}", fontsize=6, color=S.OBS, va="bottom")
    axs[0].set_yticks(y, rows.label)
    axs[0].set_ylim(-0.7, len(rows) + 0.4)
    handles = [Line2D([], [], marker="o", ls="", color=S.TRACK_COLOR["historical"], label="Concentration-driven mean"),
               Line2D([], [], marker="o", ls="", color=S.TRACK_COLOR["esm-hist"], label="Emissions-driven mean"),
               Line2D([], [], color="#999", lw=1.6, label="Member range"),
               Patch(color=S.OBS_BAND, label="GSAT ensemble 5–95%")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.55, -0.06))
    fig.tight_layout(rect=(0, 0.02, 1, 0.82))
    # Title bound to the computed consistency test (z > 1.645: one-sided 5%, obs uncertainty
    # + internal variability). Conservative toward the "too warm" claim: a model counts only
    # if every published track exceeds obs.
    t = rows.dropna(subset=["dT_2012_2021_z"]).groupby("model").dT_2012_2021_z
    n_hot, n_split = int((t.min() > 1.645).sum()), int(((t.max() > 1.645) & (t.min() <= 1.645)).sum())
    untested = sorted(set(rows.model) - set(t.groups))
    split = ", ".join(t.min()[(t.max() > 1.645) & (t.min() <= 1.645)].index)
    ts = rows.dropna(subset=["dT_2012_2021_driftcorr_seg_z"]).groupby("model").dT_2012_2021_driftcorr_seg_z
    n_hot_seg = int((ts.min() > 1.645).sum())
    S.titles(fig, f"{n_hot} of {t.ngroups} CMIP7 models tested warm more than observed over 2012–2021",
             "Ensemble mean (dot) and member range (bar) of tas; GSAT ensemble median (line) and 5–95% (shading).\n"
             "Counted if every track exceeds obs at one-sided 5% (obs uncertainty + internal variability)"
             + (f"; {split} exceeds in one track only" if n_split else "")
             + f".\nNot tested (no control run): {', '.join(untested)}. With each member's contemporaneous "
             f"control trend removed: {n_hot_seg} of {ts.ngroups}.", gap=0.06)
    S.save(fig, FIG / "fig02_historical_metrics")


def fig02b(h6):
    """CMIP6 version of fig02 from outputs/cmip6_historical_metrics.csv (scripts/cmip6_historical.py).
    One member per model and no control runs, so no z-test and no count in the title."""
    obs = h6[h6.track == "obs"].iloc[0]
    rows = h6[h6.track != "obs"].copy()
    rows["pred"] = rows.predecessor_of.fillna("") != ""
    rows["label"] = np.where(rows.pred, rows.model + "  (→ " + rows.predecessor_of.fillna("") + ")", rows.model)
    rows = rows.sort_values("dT_2012_2021")
    y = np.arange(len(rows))
    c_all, c_pred = S.SLOTS[0], S.SLOTS[2]
    fig, axs = plt.subplots(1, 2, figsize=(S.WIDTH, 175 * S.MM), sharey=True)
    for ax, k, lab in [(axs[0], "dT_2012_2021", "2012–2021 warming (°C vs 1850–1900)"),
                       (axs[1], "trend_1980_2021", "1980–2021 trend (°C per decade)")]:
        ax.axvspan(obs[f"{k}_p05"], obs[f"{k}_p95"], color=S.OBS_BAND, alpha=0.6, lw=0)
        ax.axvline(obs[k], color=S.OBS, lw=0.9)
        med = rows[k].median()
        ax.axvline(med, color=S.INK2, lw=0.8, ls=(0, (3, 2)))
        for yi, (_, r) in zip(y, rows.iterrows()):
            c = c_pred if r.pred else c_all
            ax.plot([r[f"{k}_min"], r[f"{k}_max"]], [yi, yi], color=c, lw=1.6, alpha=0.45, solid_capstyle="round")
            ax.plot(r[k], yi, "o", ms=4.5 if r.pred else 4, color=c, mec="white", mew=0.6)
        ax.set_xlabel(lab)
        ax.grid(axis="y", visible=False)
        box = dict(facecolor="white", edgecolor="none", pad=0.6)
        ax.text(obs[k], len(rows) - 0.3, f" GSAT ensemble {obs[k]:.2f}", fontsize=6, color=S.OBS, va="bottom",
                bbox=box)
        # median label on the side away from the obs line
        ax.text(med, len(rows) + 0.9, f" CMIP6 median {med:.2f} " , fontsize=6, color=S.INK2, va="bottom",
                ha="left" if med >= obs[k] else "right", bbox=box)
    axs[0].set_yticks(y, rows.label)
    for tl, p in zip(axs[0].get_yticklabels(), rows.pred):
        if p:
            tl.set_fontweight("bold")
    axs[0].set_ylim(-0.8, len(rows) + 2.2)
    handles = [Line2D([], [], marker="o", ls="", color=c_all, label="CMIP6 model (one member)"),
               Line2D([], [], marker="o", ls="", color=c_pred, label="Predecessor of a CMIP7 model"),
               Line2D([], [], color="#999", lw=1.6, label="Range across SSP splices"),
               Line2D([], [], color=S.INK2, lw=0.8, ls=(0, (3, 2)), label="CMIP6 median"),
               Patch(color=S.OBS_BAND, label="GSAT ensemble 5–95%")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.035))
    fig.tight_layout(rect=(0, 0.01, 1, 0.885))
    n_spl = rows.n_splices.astype(int)
    S.titles(fig, f"CMIP6 historical warming against the GSAT ensemble ({len(rows)} models)",
             f"One member per model from Hausfather et al. (2022), historical spliced to SSP2-4.5 from 2015 (dot); bar spans "
             f"the {n_spl.min()}–{n_spl.max()} available\nSSP splices. GSAT ensemble median (line) and 5–95% (shading). "
             "Descriptive only: with one member and no control run per model, no\ninternal-variability test is applied. "
             "HadGEM3-GC31-MM, the predecessor of UKCM2a-0-HH, is not in the dataset.", gap=0.03)
    S.save(fig, FIG / "fig02b_historical_metrics_cmip6")


def fig03(obs, pm):
    scen = [("vl", "Very low (scen7-vl)"), ("m", "Medium (scen7-m)"), ("h", "High (scen7-h)")]
    fig, axs = plt.subplots(1, 3, figsize=(S.WIDTH, 80 * S.MM), sharey=True)
    used = set()
    for ax, (s, title) in zip(axs, scen):
        ax.plot(obs.loc[:OBS_END].index, obs.loc[:OBS_END, "median"], color=S.OBS, lw=0.9)
        for m in MODELS:
            for exp, ls in ((f"scen7-{s}", "-"), (f"esm-scen7-{s}", (0, (3, 1.5)))):
                d = load(m, exp)
                if d is None:
                    continue
                c = S.MODEL_COLOR[m]
                ax.fill_between(d.index, d.min(axis=1), d.max(axis=1), color=c, alpha=0.12, lw=0)
                ax.plot(d.index, d.mean(axis=1), color=c, lw=1.1, ls=ls)
                used.add((m, ls == "-"))
        ax.set_title(title, loc="left", fontweight="bold", pad=3)
        ax.set_xlim(1950, 2100)
        ax.set_ylim(-0.5, 5.5)
        ax.axhline(0, color=S.INK2, lw=0.4)
    axs[0].set_ylabel("°C vs 1850–1900")
    mods = [m for m in MODELS if any(u[0] == m for u in used)]
    handles = [Line2D([], [], color=S.MODEL_COLOR[m], lw=1.6, label=m) for m in mods]
    handles += [Line2D([], [], color="#777", lw=1.2, label="Concentration-driven"),
                Line2D([], [], color="#777", lw=1.2, ls=(0, (3, 1.5)), label="Emissions-driven"),
                Line2D([], [], color=S.OBS, lw=1.2, label="GSAT ensemble median")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout(rect=(0, 0.02, 1, 0.84))
    S.titles(fig, "CMIP7 projections published so far",
             "Ensemble-mean tas (line) and member range (shading), °C relative to 1850–1900, by ScenarioMIP "
             "CMIP7 scenario.", gap=0.07)
    S.save(fig, FIG / "fig03_projections")


def main():
    S.setup()
    FIG.mkdir(exist_ok=True)
    obs = pd.read_csv(OUT / "obs_gsat_ensemble_summary.csv", index_col=0)
    hm = pd.read_csv(OUT / "historical_metrics.csv")
    pm = pd.read_csv(OUT / "projection_metrics.csv")
    fig01(obs, hm)
    fig02(hm)
    fig02b(pd.read_csv(OUT / "cmip6_historical_metrics.csv"))
    fig03(obs, pm)


if __name__ == "__main__":
    main()
