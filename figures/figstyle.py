"""Shared figure style: Okabe-Ito colorblind-safe palette, Optics
Express-friendly sizing (85 mm single column, 175 mm double column)."""
import pathlib
import matplotlib as mpl
from matplotlib import font_manager

# Times New Roman for all figure text (user-supplied TTFs in fonts/)
_FDIR = pathlib.Path(__file__).resolve().parent / "fonts"
for _f in ("times.ttf", "timesbd.ttf", "timesbi.ttf", "timesi.ttf"):
    _p = _FDIR / _f
    if _p.exists():
        font_manager.fontManager.addfont(str(_p))

# Okabe-Ito
OI = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9",
          green="#009E73", yellow="#F0E442", blue="#0072B2",
          vermilion="#D55E00", purple="#CC79A7", grey="#999999")

MM = 1 / 25.4
SINGLE = 85 * MM
# Full text width of the optica-article class (\linewidth = 379.42 pt),
# so text placed in figures prints at its nominal point size.
DOUBLE = 379.41753 / 72.27

def setup():
    mpl.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "font.family": "Times New Roman", "mathtext.fontset": "stix",
        "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "lines.linewidth": 1.2,
        "figure.dpi": 200, "savefig.dpi": 400,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "legend.frameon": False,
    })

def panel_label(ax, s, dx=-0.12, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="top", ha="left")
