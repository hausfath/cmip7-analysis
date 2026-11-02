"""One-off: every CMIP7 emissions-driven historical member (esm-hist), continued with esm-scen7-m
to 2026 where published, against the GSAT ensemble. Not part of the report.

  python -I scripts/oneoff_esmhist_scen7m.py  ->  figures/oneoff_esmhist_scen7m.{png,pdf}

esm-scen7-m files in derived/gmst already hold each member's full parent history (joined by parent
metadata, gmst.py), so they are used whole; esm-hist members with no esm-scen7-m child stop in 2021.
Ensemble means use only the members plotted over the whole period, so composition is fixed.
"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import figstyle as S  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
GD, OUT, FIG = ROOT / "derived" / "gmst", ROOT / "outputs", ROOT / "figures"
MODELS = ["CanESM6-0-MR", "CanESM5-1", "UKESM1-3-LL", "ACCESS-ESM1-6", "EC-Earth3-ESM-1-1"]
END, OBS_END = 2026, 2025


def members(model):
    """(full-period members to END, esm-hist-only members to 2021)."""
    h = pd.read_csv(GD / f"{model}.esm-hist.csv", index_col=0).loc[:2021]
    f = GD / f"{model}.esm-scen7-m.csv"
    if not f.exists():
        return None, h
    m = pd.read_csv(f, index_col=0).loc[:END]
    sc = pd.read_csv(OUT / "splice_check.csv")
    parents = set(sc[(sc.model == model) & (sc.scenario == "esm-scen7-m") & sc.used].parent_member)
    return m, h[[c for c in h.columns if c not in parents]]


def main():
    S.setup()
    obs = pd.read_csv(OUT / "obs_gsat_ensemble_summary.csv", index_col=0).loc[:OBS_END]
    fig, ax = plt.subplots(figsize=(S.WIDTH, 100 * S.MM))
    handles, notes, ends = [], [], []
    for model in MODELS:
        c = S.MODEL_COLOR[model]
        full, hist_only = members(model)
        main_ = full if full is not None else hist_only
        for d, z in ((main_, 2), (hist_only if full is not None else None, 2)):
            if d is not None and d.shape[1]:
                ax.plot(d.index, d.to_numpy(), color=c, lw=0.35, alpha=0.35, zorder=z)
        mean = main_.mean(axis=1).dropna()
        ax.plot(mean.index, mean, color=c, lw=1.3, zorder=4)
        ends.append([mean.iloc[-1], model, c, mean.index[-1]])
        n_extra = 0 if full is None else hist_only.shape[1]
        lab = f"{model} (n={main_.shape[1]}{', to 2021' if full is None else ''})"
        handles.append(Line2D([], [], color=c, lw=1.6, label=lab))
        if n_extra:
            notes.append(f"{n_extra} {model} members without an esm-scen7-m run are drawn to 2021 and left out of its mean")
    ax.fill_between(obs.index, obs.p05, obs.p95, color=S.OBS_BAND, alpha=0.7, lw=0, zorder=3)
    ax.plot(obs.index, obs["median"], color=S.OBS, lw=1.2, zorder=5)
    ends.append([obs["median"].iloc[-1], f"GSAT ensemble ({OBS_END})", S.OBS, OBS_END])
    # direct labels at the right edge, nudged apart (min gap in °C), joined to each line end
    ends.sort()
    GAP = 0.14
    ys = [e[0] for e in ends]
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + GAP)
    for (yv, name, c, x_end), yl in zip(ends, ys):
        ax.plot([x_end + 0.3, END + 1.6], [yv, yl], color=c, lw=0.5, zorder=4)
        ax.text(END + 1.9, yl, name, fontsize=6, color=S.OBS, va="center", ha="left")
    ax.axhline(0, color=S.INK2, lw=0.4)
    ax.axvline(2021.5, color=S.INK2, lw=0.5, ls=(0, (2, 2)))
    ax.text(2021.3, ax.get_ylim()[0] + 0.05, "esm-hist | esm-scen7-m ", fontsize=6, color=S.INK2, ha="right", va="bottom")
    ax.set_xlim(1850, END + 22)
    ax.set_xticks(range(1860, END + 1, 20))
    ax.spines["bottom"].set_bounds(1850, END + 1)
    ax.set_ylabel("°C vs 1850–1900")
    handles += [Line2D([], [], color=S.OBS, lw=1.2, label=f"GSAT ensemble median (to {OBS_END})"),
                Patch(color=S.OBS_BAND, label="GSAT ensemble 5–95%"),
                Line2D([], [], color="#999", lw=0.6, label="Individual members")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.11))
    fig.tight_layout(rect=(0, 0.02, 1, 0.83))
    sub = ("Annual tas, each member relative to its own 1850–1900 mean; thick lines are ensemble means. "
           "Continued with the medium\nscenario (esm-scen7-m) from 2022 where published; ACCESS-ESM1-6 and EC-Earth3-ESM-1-1 "
           "have no esm-scen7-m run yet.")
    if notes:
        sub += "\n" + "; ".join(notes) + "."
    S.titles(fig, f"CMIP7 emissions-driven runs against observations, 1850–{END}", sub, gap=0.065)
    S.save(fig, FIG / "oneoff_esmhist_scen7m")


if __name__ == "__main__":
    main()
