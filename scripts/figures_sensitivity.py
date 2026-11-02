"""Climate-sensitivity figures from outputs/*_sensitivity.csv and derived/gregory (no raw data).

  python -I scripts/figures_sensitivity.py

fig04_ecs_tcr_vs_cmip6   ECS and TCR of each CMIP7 model on the CMIP6 distribution, with AR6
                         assessed ranges and an arrow from the nearest CMIP6 predecessor
fig05_gregory            Gregory plots, CMIP7 vs predecessor
fig06_ecs_tcr_scatter    ECS vs TCR
Also writes outputs/lineage_pairs.csv.
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
OUT, FIG, GDIR = ROOT / "outputs", ROOT / "figures", ROOT / "derived" / "gregory"
# Nearest CMIP6 predecessor by model family (our assumption; not an official lineage table).
LINEAGE = {"CanESM6-0-MR": ("CanESM5", "r1i1p1f1"), "CanESM5-1": ("CanESM5", "r1i1p2f1"),
           "UKESM1-3-LL": ("UKESM1-0-LL", "r1i1p1f2"), "UKCM2-0-LL": ("HadGEM3-GC31-LL", "r1i1p1f3"),
           "UKCM2a-0-HH": ("HadGEM3-GC31-MM", "r1i1p1f3"), "ACCESS-ESM1-6": ("ACCESS-ESM1-5", "r1i1p1f1")}
AR6 = {"ECS": (3.0, (2.5, 4.0), (2.0, 5.0)), "TCR": (1.8, (1.4, 2.2), (1.2, 2.4))}
C7 = S.SLOTS[0]
C6 = "#8d8c86"
LABEL_OFFSET = {"UKESM1-3-LL": (0, -12, "center"), "CanESM6-0-MR": (0, 7, "center"), "CanESM5-1": (6, -8, "left")}


def pairs():
    c7 = pd.read_csv(OUT / "cmip7_sensitivity.csv")
    c6 = pd.read_csv(OUT / "cmip6_sensitivity.csv")
    ref = pd.read_csv(ROOT / "reference" / "hausfather2022_cmip6_ecs_tcr.csv", comment="#")
    ref = ref[~ref.copy_of_other_model].set_index("model")
    rows = []
    for _, r in c7.iterrows():
        pm, pmem = LINEAGE[r.model]
        p = c6[(c6.model == pm) & (c6.member == pmem)]
        p = p.iloc[0] if len(p) else None
        rows.append({"model": r.model, "member": r.member, "ECS": r.ECS_rm21, "ECS_p05": r.ECS_p05_rm21,
                     "ECS_p95": r.ECS_p95_rm21, "TCR": r.get("TCR_linear", np.nan),
                     "F2x": r.F2x_rm21, "lambda": r.lambda_rm21, "ECS_21_150": r.ECS_21_150_rm21,
                     "predecessor": pm, "pred_member": pmem,
                     # same-code value where the predecessor was recomputed, else the published one
                     "pred_ECS": p.ECS_rm21 if p is not None else ref.loc[pm, "ECS150"],
                     "pred_ECS_p05": p.ECS_p05_rm21 if p is not None else np.nan,
                     "pred_ECS_p95": p.ECS_p95_rm21 if p is not None else np.nan,
                     "pred_TCR": p.TCR_linear if p is not None else ref.loc[pm, "TCR"],
                     "pred_F2x": p.F2x_rm21 if p is not None else np.nan,
                     "pred_lambda": p.lambda_rm21 if p is not None else np.nan,
                     "pred_source": "recomputed" if p is not None else "H22 published"})
    df = pd.DataFrame(rows)
    df["dECS"] = df.ECS - df.pred_ECS
    df["dTCR"] = df.TCR - df.pred_TCR
    # Headline (conservative): "lower/higher" only when the two 5–95 % regression ranges do not
    # overlap. Sensitivity: predecessor outside the CMIP7 range only. Neither range includes
    # realisation-to-realisation ECS spread. A missing predecessor is "untested", never "no change".
    def change(lo_new, hi_new, lo_old, hi_old):
        out = np.where(lo_old > hi_new, "lower", np.where(hi_old < lo_new, "higher", "no clear change"))
        return np.where(np.isfinite(lo_old) & np.isfinite(hi_old), out, "untested")
    df["ECS_change"] = change(df.ECS_p05, df.ECS_p95, df.pred_ECS_p05, df.pred_ECS_p95)
    df["ECS_change_single_range"] = change(df.ECS_p05, df.ECS_p95, df.pred_ECS, df.pred_ECS)
    df.to_csv(OUT / "lineage_pairs.csv", index=False)
    return df, ref


def ar6_bands(ax, key, label=True):
    best, lik, vl = AR6[key]
    ax.axvspan(*vl, color="#ecebe6", lw=0, zorder=0)
    ax.axvspan(*lik, color="#d9d8d1", lw=0, zorder=0)
    ax.axvline(best, color=S.INK2, lw=0.6, ls=(0, (2, 2)), zorder=1)


def fig04(df, ref):
    df = df.sort_values("ECS")
    fig, axs = plt.subplots(1, 2, figsize=(S.WIDTH, 92 * S.MM), sharey=True)
    n6 = {"ECS": int(ref.ECS150.notna().sum()), "TCR": int(ref.TCR.notna().sum())}
    rng = np.random.default_rng(1)
    for ax, key, col6 in [(axs[0], "ECS", "ECS150"), (axs[1], "TCR", "TCR")]:
        ar6_bands(ax, key)
        v6 = ref[col6].dropna().to_numpy()
        top = len(df) + 0.9
        ax.scatter(v6, top + rng.uniform(-0.32, 0.32, len(v6)), s=9, color=C6, alpha=0.75, lw=0, zorder=2)
        ax.plot(np.median(v6), top, "|", ms=14, mew=1.4, color=S.OBS, zorder=3)
        for y, (_, r) in enumerate(df.iterrows()):
            pv, v = r[f"pred_{key}"], r[key]
            if not np.isfinite(v):
                ax.text(pv + 0.06, y, "no CMIP7 1pctCO2 run yet", fontsize=6, va="center", color=S.INK2)
            if np.isfinite(pv) and abs(v - pv) > 0.15:   # shorter arrows vanish under the markers
              ax.annotate("", xy=(v, y), xytext=(pv, y), zorder=3,
                        arrowprops=dict(arrowstyle="-|>", color=C6, lw=0.8, shrinkA=3.5, shrinkB=4.5,
                                        mutation_scale=7))
            ax.plot(pv, y, "o", ms=5, mfc="white", mec=C6, mew=1.0, zorder=4)
            if key == "ECS":
                ax.plot([r.ECS_p05, r.ECS_p95], [y, y], color=C7, lw=2.2, alpha=0.5, solid_capstyle="round",
                        zorder=4)
            ax.plot(v, y, "o", ms=6, color=C7, mec="white", mew=0.6, zorder=5)
            ax.text(v, y + 0.3, f"{v:.2f}", fontsize=6, ha="center", va="bottom", color=S.OBS, zorder=6)
        ax.set_xlabel(f"{key} (°C)")
        ax.grid(axis="y", visible=False)
    ylabels = [f"{r.model}\n← {r.predecessor}" for _, r in df.iterrows()] + [f"CMIP6 (n={n6['ECS']}/{n6['TCR']})"]
    axs[0].set_yticks(list(range(len(df))) + [len(df) + 0.9], ylabels, fontsize=6.5)
    axs[0].set_ylim(-0.7, len(df) + 1.5)
    axs[0].set_xlim(1.5, 6.2)
    axs[1].set_xlim(1.0, 3.1)
    nb = int((df.ECS_change == "lower").sum())
    nb1 = int((df.ECS_change_single_range == "lower").sum())
    nun = int((df.ECS_change == "untested").sum())
    above_lik = int((df.ECS > AR6["ECS"][1][1]).sum())
    above_vl = int((df.ECS > AR6["ECS"][2][1]).sum())
    handles = [Line2D([], [], marker="o", ls="", color=C7, mec="white", label="CMIP7 (ECS bar: 5–95% regression range)"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=C6, label="Nearest CMIP6 predecessor (same code)"),
               Line2D([], [], marker="o", ls="", color=C6, ms=4, label="CMIP6 models (Hausfather et al. 2022)"),
               Line2D([], [], marker="|", ls="", color=S.OBS, ms=9, mew=1.4, label="CMIP6 median"),
               Patch(color="#d9d8d1", label="AR6 likely"), Patch(color="#ecebe6", label="AR6 very likely"),
               Line2D([], [], color=S.INK2, lw=0.8, ls=(0, (2, 2)), label="AR6 best estimate")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout(rect=(0, 0.03, 1, 0.80))
    S.titles(fig, f"{nb} of {len(df)} CMIP7 models so far have lower ECS than their CMIP6 predecessor",
             f"But {above_lik} of {len(df)} remain above the AR6 likely ECS range (2.5–4 °C) and {above_vl} above the "
             f"very likely range (2–5 °C).\nECS: Gregory regression, 150 yr abrupt-4xCO2. TCR: 1pctCO2 years 61–80. "
             f"Lower = 5–95% regression ranges do not overlap ({nb1} if only the\npredecessor's central value must "
             f"lie outside the CMIP7 range" + (f"; {nun} pair untested" if nun else "") + "). Predecessor pairings by "
             "model family (our assumption).", gap=0.06)
    S.save(fig, FIG / "fig04_ecs_tcr_vs_cmip6")


def fig05(df):
    df = df.sort_values("ECS", ascending=False)
    fig, axs = plt.subplots(2, 3, figsize=(S.WIDTH, 118 * S.MM), sharex=True, sharey=True)
    for ax, (_, r) in zip(axs.flat, df.iterrows()):
        for model, mem, c, mfc, lab in [(r.predecessor, r.pred_member, C6, "none", r.predecessor),
                                        (r.model, r.member, C7, C7, r.model)]:
            f = GDIR / f"{model}_{mem}.csv"
            if not f.exists():
                continue
            g = pd.read_csv(f)
            ax.scatter(g.dT, g.dN, s=5, facecolors=mfc, edgecolors=c, lw=0.5, alpha=0.8, zorder=3)
            b, a = np.polyfit(g.dT, g.dN, 1)
            x = np.array([0, -a / b])
            ax.plot(x, a + b * x, color=c, lw=1.1, zorder=4)
        ax.axhline(0, color=S.INK2, lw=0.5)
        ax.set_title(r.model, loc="left", fontweight="bold", pad=3)
        ax.text(0.97, 0.95, f"ECS {r.ECS:.2f} °C\n{r.predecessor}: " + (f"{r.pred_ECS:.2f} °C" if np.isfinite(r.pred_ECS) else "n/a"), transform=ax.transAxes,
                ha="right", va="top", fontsize=6.3, color=S.INK2, linespacing=1.4)
        ax.set_xlim(0, 12.5)
        ax.set_ylim(-1.5, 9.5)
    for ax in axs[:, 0]:
        ax.set_ylabel("ΔN (W m$^{-2}$)")
    for ax in axs[-1]:
        ax.set_xlabel("ΔT (°C)")
    handles = [Line2D([], [], marker="o", ls="-", color=C7, ms=3, label="CMIP7 model"),
               Line2D([], [], marker="o", ls="-", color=C6, mfc="none", ms=3, label="Nearest CMIP6 predecessor")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.02, 1, 0.88))
    S.titles(fig, "Gregory regressions for CMIP7 abrupt-4xCO2 runs",
             "Annual-mean TOA imbalance vs surface warming, years 1–150, anomalies vs a 21-yr running mean of "
             "the contemporaneous piControl.\nThe x-intercept of each line is twice the ECS.", gap=0.05)
    S.save(fig, FIG / "fig05_gregory")


def fig06(df, ref):
    fig, ax = plt.subplots(figsize=(S.WIDTH * 0.62, 95 * S.MM))
    ar = AR6
    ax.add_patch(plt.Rectangle((ar["ECS"][1][0], ar["TCR"][1][0]), np.diff(ar["ECS"][1])[0],
                               np.diff(ar["TCR"][1])[0], color="#d9d8d1", lw=0, zorder=0))
    r6 = ref.dropna(subset=["ECS150", "TCR"])
    ax.scatter(r6.ECS150, r6.TCR, s=10, color=C6, alpha=0.7, lw=0, zorder=2)
    for _, r in df.dropna(subset=["TCR"]).iterrows():
        ax.plot(r.ECS, r.TCR, "o", ms=6, color=C7, mec="white", mew=0.6, zorder=4)
        dx, dy, ha = LABEL_OFFSET.get(r.model, (5, -2, "left"))
        ax.annotate(r.model, (r.ECS, r.TCR), xytext=(dx, dy), textcoords="offset points", fontsize=6.3,
                    zorder=5, ha=ha)
    ax.set_xlabel("ECS (°C)")
    ax.set_ylabel("TCR (°C)")
    ax.set_xlim(1.5, 6.2)
    ax.set_ylim(1.0, 3.1)
    handles = [Line2D([], [], marker="o", ls="", color=C7, label="CMIP7"),
               Line2D([], [], marker="o", ls="", color=C6, ms=4, label="CMIP6 models (Hausfather et al. 2022)"),
               Patch(color="#d9d8d1", label="AR6 likely ECS × TCR")]
    ax.legend(handles=handles, loc="upper left", frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    S.titles(fig, "ECS and TCR of CMIP7 models so far",
             "CMIP7 models placed on the CMIP6 distribution (predecessor changes: Fig. 4).\nUKCM2a-0-HH has no 1pctCO2 run yet.", gap=0.075)
    S.save(fig, FIG / "fig06_ecs_tcr_scatter")


def main():
    S.setup()
    df, ref = pairs()
    with pd.option_context("display.width", 250):
        print(df[["model", "ECS", "ECS_p05", "ECS_p95", "pred_ECS", "dECS", "ECS_change", "TCR", "pred_TCR", "dTCR",
                  "pred_source"]].round(2).to_string())
    fig04(df, ref)
    fig05(df)
    fig06(df, ref)


if __name__ == "__main__":
    main()
