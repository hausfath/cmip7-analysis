"""Figures for the Earth energy imbalance comparison (scripts/eei.py).

fig07_eei_timeseries   one panel per model: 10-yr centred rates, models vs observations
fig08_eei_windows      ensemble means and members for 2005-2019, 2015-2024 and the change
"""
import pathlib, sys
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import figstyle as S  # noqa: E402
import eei as E  # noqa: E402

plt = S.plt
FIG = E.L.ROOT / "figures"
WIN = 10                      # running window length (yr)
Y0, Y1 = 1960, 2024
OBS_STYLE = {"wu": ("Wu et al. 2025 (ocean heat)", S.SLOTS[2]), "miniere": ("Minière et al. 2026 (ocean heat)", S.SLOTS[6]),
             "ceres": ("CERES EBAF (anchored to in-situ 2005–2015)", S.OBS)}
TRACK_LS = {"historical": "-", "esm-hist": (0, (3, 1.5))}
OBS_DY = 1.4                  # vertical spacing of observational rows in fig08
TRACK_SHORT = {"historical": "conc.", "esm-hist": "emis."}


def running(fn, years):
    return pd.Series({y0 + (WIN - 1) / 2: fn(y0, y0 + WIN - 1) for y0 in years})


def obs_running():
    wu = pd.read_csv(E.OBS / "wu2025_gohc_0-2000m_ZJ.csv", index_col=0) * 1e21 / E.A_EARTH
    mi = pd.read_csv(E.OBS / "miniere2026_gohc_0-2000m_J_no_mask.csv", index_col=0).drop(columns="ensemble_mean") / E.A_EARTH
    ce = pd.read_csv(E.OBS / "ceres_ebaf_ed4.2.1_annual.csv", index_col=0).net_toa_Wm2
    out = {}
    for name, df in [("wu", wu), ("miniere", mi)]:
        rows = {}
        for y0 in range(Y0, Y1 - WIN + 2):
            r = df.apply(lambda c: E.ols_rate(c, y0, y0 + WIN - 1)).dropna()
            if len(r) >= (10 if name == "wu" else 3):
                e = (r + E.deep(y0, y0 + WIN - 1)) / E.F_OCEAN
                rows[y0 + (WIN - 1) / 2] = (e.mean(), e.quantile(0.05) if name == "wu" else e.min(),
                                           e.quantile(0.95) if name == "wu" else e.max())
        out[name] = pd.DataFrame(rows, index=["mid", "lo", "hi"]).T
    out["ceres"] = running(lambda a, b: ce.loc[a:b].mean() if ce.loc[a:b].size == WIN else np.nan,
                           range(ce.index.min(), Y1 - WIN + 2)).dropna()
    return out


def model_running(df):
    """10-yr OLS rate of cumulative heat per member (the estimator applied to OHC), then the mean."""
    res = {}
    for m in df:
        n = df[m].dropna()
        h = E.cum_heat(n)
        res[m] = running(lambda a, b: E.ols_rate(h, a, b), range(Y0, int(n.index.max()) - WIN + 2))
    r = pd.DataFrame(res)
    return r.mean(axis=1), r.min(axis=1), r.max(axis=1)


def fig07(met):
    obs = obs_running()
    mods = sorted(met.model.unique(), key=lambda m: -met[met.model == m].rate_2005_2019.mean())
    nc = 3
    nr = -(-len(mods) // nc)
    fig, axs = plt.subplots(nr, nc, figsize=(S.WIDTH, (30 + 42 * nr) * S.MM), sharex=True, sharey=True, squeeze=False)
    for ax, model in zip(axs.flat, mods):
        o = obs["wu"]
        ax.fill_between(o.index, o.lo, o.hi, color=OBS_STYLE["wu"][1], alpha=0.18, lw=0)
        ax.plot(o.index, o.mid, color=OBS_STYLE["wu"][1], lw=0.9)
        o = obs["miniere"]
        ax.fill_between(o.index, o.lo, o.hi, color=OBS_STYLE["miniere"][1], alpha=0.18, lw=0)
        ax.plot(o.index, o.mid, color=OBS_STYLE["miniere"][1], lw=0.9)
        ax.plot(obs["ceres"].index, obs["ceres"].values, color=S.OBS, lw=1.1, ls=(0, (1, 1)))
        for _, r in met[met.model == model].iterrows():
            df = pd.read_csv(E.DER / f"{model}.{r.track}.csv", index_col=0)
            mid, lo, hi = model_running(df)
            c = S.TRACK_COLOR[r.track]
            if df.shape[1] > 1:
                ax.fill_between(mid.index, lo, hi, color=c, alpha=0.15, lw=0)
            ax.plot(mid.index, mid, color=c, lw=1.3, ls=TRACK_LS[r.track])
        n = ", ".join(f"{'conc.' if r.track == 'historical' else 'emis.'} n={r.n_members}"
                      for _, r in met[met.model == model].iterrows())
        ax.set_title(f"{model}\n", fontsize=7.5, loc="left", fontweight="bold")
        ax.text(0, 1.02, n, transform=ax.transAxes, fontsize=6, color=S.INK2, va="bottom")
        ax.axhline(0, color=S.INK2, lw=0.5)
        ax.set_xlim(Y0 + WIN / 2 - 1, Y1 - WIN / 2 + 1)
    for ax in axs[:, 0]:
        ax.set_ylabel("EEI (W m$^{-2}$)")
    for ax in axs.flat[len(mods):]:
        ax.set_visible(False)
    for ax in axs[-1]:
        ax.set_xlabel(f"Centre of {WIN}-yr window")
    for ax in axs[-2, len(mods) - nc * (nr - 1):] if nr > 1 and len(mods) % nc else []:
        ax.tick_params(labelbottom=True)
        ax.set_xlabel(f"Centre of {WIN}-yr window")
    handles = [Line2D([], [], color=S.TRACK_COLOR["historical"], lw=1.3, label="Concentration-driven, ensemble mean"),
               Line2D([], [], color=S.TRACK_COLOR["esm-hist"], lw=1.3, ls=TRACK_LS["esm-hist"], label="Emissions-driven, ensemble mean"),
               Patch(color="#999", alpha=0.35, label="Member range"),
               Line2D([], [], color=OBS_STYLE["wu"][1], lw=0.9, label=OBS_STYLE["wu"][0] + ", 5–95%"),
               Line2D([], [], color=OBS_STYLE["miniere"][1], lw=0.9, label=OBS_STYLE["miniere"][0] + ", product range"),
               Line2D([], [], color=S.OBS, lw=1.1, ls=(0, (1, 1)), label=OBS_STYLE["ceres"][0])]
    spare = list(axs.flat[len(mods):])
    if spare:
        spare[0].set_visible(True)
        spare[0].axis("off")
        spare[0].legend(handles=handles, loc="center left", frameon=False, fontsize=6)
    else:
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0, 1, 1 - 23 / (30 + 42 * nr)))
    S.titles(fig, "Earth's energy imbalance in CMIP7 models and observations",
             f"Rate of heat gain over centred {WIN}-yr windows (W m$^{{-2}}$ of Earth's surface): OLS slope of cumulative heat for "
             "models and the ocean-heat\nestimates, window mean for CERES. Ocean heat 0–2000 m plus the deep ocean "
             f"({E.DEEP} W m$^{{-2}}$ from {E.DEEP_FROM}), divided by the ocean's {E.F_OCEAN:.0%} share of\nEarth's heat gain. "
             "Model net TOA flux with the control run's linear fit removed; windows after 2021 use only the members joined to a scenario run.", gap=0.05 * 160 / (30 + 42 * nr))
    S.save(fig, FIG / "fig07_eei_timeseries")


def fig08(met, obs, tests):
    panels = [("rate_2005_2019", "2005–2019 (W m$^{-2}$)"), ("rate_2015_2024", "2015–2024 (W m$^{-2}$)"),
              ("change", "2015–2024 minus 2005–2014 (W m$^{-2}$)")]
    rows = met.copy()
    rows["label"] = rows.model + "  (" + rows.track.map(TRACK_SHORT) + \
        ", n=" + rows.n_members.astype(str) + ")"
    rows = rows.sort_values("rate_2005_2019")
    names = {"wu": "Wu 2025", "miniere": "Minière 2026", "hakuba_hydrography": "Hakuba 2021 hydr.",
             "hakuba_geodetic": "Hakuba 2021 geod.", "ceres": "CERES EBAF"}
    fig, axs = plt.subplots(1, 3, figsize=(S.WIDTH, 118 * S.MM), sharey=True)
    nmod = len(rows)
    for ax, (k, lab) in zip(axs, panels):
        o = obs[obs.metric == k]
        if k.startswith("rate_"):   # CERES level shown as window mean, for reference
            o = pd.concat([o, obs[obs.metric == k.replace("rate_", "mean_")]])
        for j, r in enumerate(o.itertuples()):
            yy = nmod + 0.6 + OBS_DY * j
            c = S.OBS
            if pd.notna(getattr(r, "p05", np.nan)):
                ax.plot([r.p05, r.p95], [yy, yy], color=c, lw=1.4, alpha=0.5, solid_capstyle="round")
            ax.plot(r.value, yy, "D", ms=3.6, color=c if r.source != "ceres" or k == "change" else "white",
                    mec=c, mew=0.8)
        for i, (_, r) in enumerate(rows.iterrows()):
            c = S.TRACK_COLOR[r.track]
            df = pd.read_csv(E.DER / f"{r.model}.{r.track}.csv", index_col=0)
            vals = [E.model_metrics(df[m])[k] for m in df]
            vals = [v for v in vals if pd.notna(v)]
            ax.plot(vals, [i] * len(vals), "o", ms=2.2, color=c, alpha=0.45, mew=0)
            if vals:
                ax.plot(np.mean(vals), i, "o", ms=4.8, color=c, mec="white", mew=0.6)
                if len(vals) != r.n_members:   # fewer members reach this window (scenario joins)
                    ax.annotate(f"n={len(vals)}", (max(vals), i), xytext=(4, 0), textcoords="offset points",
                                fontsize=6, va="center", color=S.INK2)
        ax.axhline(nmod - 0.2, color=S.GRID, lw=0.8)
        ax.set_xlabel(lab)
        ax.grid(axis="y", visible=False)
    labels = list(rows.label)
    ymax = 0
    for ax, (k, _) in zip(axs, panels):
        o = obs[obs.metric == k]
        if k.startswith("rate_"):
            o = pd.concat([o, obs[obs.metric == k.replace("rate_", "mean_")]])
        ymax = max(ymax, len(o))
    axs[0].set_yticks(list(range(nmod)), labels)
    for ax, (k, _) in zip(axs, panels):
        o = obs[obs.metric == k]
        if k.startswith("rate_"):
            o = pd.concat([o, obs[obs.metric == k.replace("rate_", "mean_")]])
        for j, r in enumerate(o.itertuples()):
            ax.annotate(names[r.source], (r.value, nmod + 0.6 + OBS_DY * j), xytext=(0, 3.2), textcoords="offset points",
                        fontsize=6, ha="center", va="bottom", color=S.INK2)
    axs[0].set_ylim(-0.7, nmod + 0.6 + OBS_DY * (ymax - 1) + 1.0)
    handles = [Line2D([], [], marker="o", ls="", color=S.TRACK_COLOR["historical"], label="Concentration-driven mean"),
               Line2D([], [], marker="o", ls="", color=S.TRACK_COLOR["esm-hist"], label="Emissions-driven mean"),
               Line2D([], [], marker="o", ls="", ms=2.5, color="#888", label="Members"),
               Line2D([], [], marker="D", ls="", color=S.OBS, label="Observational estimate (5–95% or range)"),
               Line2D([], [], marker="D", ls="", color="white", mec=S.OBS, label="CERES level (not tested)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout(rect=(0, 0.02, 1, 0.8))
    v = tests.drop_duplicates(["model", "track", "metric"])
    flag = v[v.verdict != "consistent with at least one"]
    n_mod = v[v.metric == "rate_2005_2019"].model.nunique()
    n_ok = n_mod - flag[flag.metric == "rate_2005_2019"].model.nunique()
    def flagged(metric):
        f = flag[flag.metric == metric]
        return "none" if f.empty else ", ".join(
            f"{r.model} {TRACK_SHORT[r.track]} ({r.verdict.split()[0]})" for r in f.itertuples())
    ce = tests[(tests.metric == "change") & (tests.obs == "ceres")]
    ce_lo, ce_n = int((ce.z < -E.Z).sum()), len(ce)
    S.titles(fig, f"{n_ok} of {n_mod} CMIP7 models match at least one observed energy imbalance for 2005–2019",
             "Ensemble mean (large dot) and members (small dots); observational estimates above (diamonds, 5–95% or product "
             "range). Rates are OLS slopes\nof cumulative heat. CERES levels are window means tied to an in-situ value, so "
             "they are shown but not tested; 2015–2024 is tested against Minière\nonly. An ensemble is flagged only if it differs "
             "from every tested estimate (|z| > 1.645, observational plus internal-variability uncertainty).\n"
             f"Flagged: 2005–2019 {flagged('rate_2005_2019')}; 2015–2024 {flagged('rate_2015_2024')}; "
             f"change {flagged('change')}. Against CERES alone, {ce_lo} of {ce_n} ensembles gain less (no CERES uncertainty included).", gap=0.055)
    S.save(fig, FIG / "fig08_eei_windows")


def main():
    S.setup()
    met = pd.read_csv(E.OUT / "eei_metrics.csv")
    obs = pd.read_csv(E.OUT / "eei_obs.csv")
    tests = pd.read_csv(E.OUT / "eei_tests.csv")
    fig07(met)
    fig08(met, obs, tests)


if __name__ == "__main__":
    main()
