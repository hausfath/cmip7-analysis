"""Shared figure style: 180 mm wide, Arial, text >= 6 pt at print size."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MM = 1 / 25.4
WIDTH = 180 * MM

# dataviz reference palette, fixed slot order (never cycled)
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MODEL_COLOR = dict(zip(["CanESM6-0-MR", "CanESM5-1", "UKESM1-3-LL", "UKCM2-0-LL", "UKCM2a-0-HH",
                        "ACCESS-ESM1-6", "CAS-ESM2-1-CSM", "EC-Earth3-ESM-1-1"], SLOTS))
TRACK_COLOR = {"historical": SLOTS[0], "esm-hist": SLOTS[1]}
TRACK_LABEL = {"historical": "Concentration-driven", "esm-hist": "Emissions-driven"}
OBS = "#0b0b0b"
OBS_BAND = "#b9b8b2"
INK2 = "#52514e"
GRID = "#e4e3df"


def setup():
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
        "axes.linewidth": 0.6, "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.axisbelow": True, "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42,
        "svg.fonttype": "none"})


def titles(fig, title, subtitle, top=0.985, gap=0.045):
    """Title and subtitle with at least one line-height between them."""
    fig.text(0.01, top, title, fontsize=10, fontweight="bold", va="top", ha="left")
    fig.text(0.01, top - gap, subtitle, fontsize=7.5, color=INK2, va="top", ha="left")


def save(fig, path):
    fig.savefig(path.with_suffix(".png"), bbox_inches="tight", facecolor="white")
    # no CreationDate, so reruns give byte-identical PDFs
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", metadata={"CreationDate": None})
    plt.close(fig)
