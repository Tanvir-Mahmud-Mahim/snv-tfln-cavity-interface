"""Fig. 1: device concept (2D top + side views with a magnified cavity
region) + level alignment + headline result."""
import json, pathlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, Ellipse

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figstyle import OI, DOUBLE, setup, panel_label
setup()
RES = HERE.parent / "sim" / "results"
frl = json.load(open(RES / "fr_vs_lambda.json"))

DIA = "#c9d4e8"      # diamond
LN = "#9ed4f0"       # TFLN (light blue so labels stay readable)
SUB = "#f0f0f0"      # SiO2
DET = "#3f3f3f"      # detector


def _lead(ax, xytext, xy, s, fs=7, color="k", ha="center", va="center",
          leader=True):
    """Upright label with an optional thin leader line to its feature."""
    arr = dict(arrowstyle="-", lw=0.5, color="#999999",
               shrinkA=2, shrinkB=1) if leader else None
    ax.annotate(s, xy=xy, xytext=xytext, fontsize=fs, color=color,
                ha=ha, va=va, zorder=30, annotation_clip=False,
                arrowprops=arr)


# hole-pattern helper: schematic pitch sequence for the overview
def overview_holes(x0=0.70):
    gaps = ([0.30]*6 + [0.290, 0.280, 0.270, 0.262, 0.258]
            + [0.258, 0.262, 0.270, 0.280, 0.290] + [0.30]*2)
    xs = [x0]
    for g in gaps:
        xs.append(xs[-1] + g)
    rr = [0.30*g for g in [0.30] + gaps]
    x_cav = x0 + 6*0.30 + sum([0.290, 0.280, 0.270, 0.262, 0.258]) + 0.258/2
    return xs, rr, x_cav


def panel_top(ax):
    """Plan (top) view of the whole device."""
    ax.set_xlim(0, 11); ax.set_ylim(0, 3.0)
    ax.set_aspect("equal"); ax.axis("off")
    yc = 1.42
    ax.add_patch(Rectangle((0, 0), 11, 3.0, fc=SUB, ec="#bbbbbb", lw=0.5))
    # TFLN bus (600 nm) ending in the detector
    ax.add_patch(Rectangle((6.05, yc-0.22), 4.25, 0.44, fc=LN, ec="k", lw=0.5))
    ax.add_patch(Rectangle((10.3, yc-0.36), 0.55, 0.72, fc=DET, ec="k",
                           lw=0.5))
    # diamond nanobeam (280 nm) + adiabatic taper onto the bus
    ax.add_patch(Rectangle((0.4, yc-0.103), 5.8, 0.206, fc=DIA, ec="k",
                           lw=0.5, zorder=5))
    ax.add_patch(Polygon([(6.2, yc-0.103), (8.0, yc), (6.2, yc+0.103)],
                         closed=True, fc=DIA, ec="k", lw=0.5, zorder=6))
    # holes
    xs, rr, x_cav = overview_holes()
    for xh, r in zip(xs, rr):
        ax.add_patch(plt.Circle((xh, yc), r, fc="white", ec="k",
                                lw=0.3, zorder=7))
    # emitter
    ax.plot([x_cav], [yc], marker="*", ms=7, color=OI["vermilion"],
            markeredgecolor="k", markeredgewidth=0.4, zorder=8)
    # photon wavepackets + propagation arrow on the bus
    for xw in (8.35, 9.1):
        tt = np.linspace(0, 1, 60)
        ax.plot(xw + 0.5*tt, yc + 0.09*np.sin(2*np.pi*3.2*tt)
                * np.exp(-((tt-0.5)/0.3)**2), color=OI["vermilion"],
                lw=1.0, zorder=8)
    ax.annotate("", xy=(10.15, yc), xytext=(9.7, yc),
                arrowprops=dict(arrowstyle="-|>", lw=1.2,
                                color=OI["vermilion"]), zorder=8)
    # static field B
    ax.annotate("", xy=(1.35, 1.02), xytext=(0.65, 0.58),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color="k"))
    ax.text(1.48, 1.02, "$B$", fontsize=8, va="center")
    # region of interest (magnified below)
    ax.add_patch(Rectangle((2.3, yc-0.33), 4.1, 0.66, fill=False,
                           ec="#777777", lw=0.7, ls=(0, (3, 2)), zorder=9))
    # labels
    _lead(ax, (1.9, 2.55), (1.9, yc+0.14), "diamond nanobeam")
    _lead(ax, (7.1, 2.55), (6.9, yc+0.07), "adiabatic taper")
    _lead(ax, (8.6, 0.68), (8.6, yc-0.24), "TFLN bus, 600 nm")
    _lead(ax, (10.35, 2.55), (10.55, yc+0.38), "detector")
    ax.text(0.22, 0.22, "SiO$_2$", fontsize=7.5, color="#666666")
    ax.text(10.85, 0.10, "top view, not to scale", fontsize=6,
            color="#999999", ha="right")


def panel_side(ax):
    """Side view along the beam: the layer stack."""
    ax.set_xlim(0, 11); ax.set_ylim(-0.02, 1.35)
    ax.set_aspect("equal"); ax.axis("off")
    # SiO2 substrate
    ax.add_patch(Rectangle((0, 0), 11, 0.46, fc=SUB, ec="k", lw=0.5))
    # TFLN film (fully etched ridge, 190 nm)
    ax.add_patch(Rectangle((6.05, 0.46), 4.25, 0.19, fc=LN, ec="k", lw=0.5))
    # diamond beam (200 nm): on SiO2, ramping onto the LN ridge
    beam = Polygon([(0.4, 0.46), (6.05, 0.46), (6.4, 0.65), (8.0, 0.65),
                    (8.0, 0.85), (6.4, 0.85), (6.05, 0.66), (0.4, 0.66)],
                   closed=True, fc=DIA, ec="k", lw=0.5, zorder=5)
    ax.add_patch(beam)
    # holes appear as through-thickness slots in side view
    xs, rr, x_cav = overview_holes()
    for xh, r in zip(xs, rr):
        ax.add_patch(Rectangle((xh-r, 0.466), 2*r, 0.188, fc="white",
                               ec="k", lw=0.25, zorder=6))
    ax.plot([x_cav], [0.56], marker="*", ms=6, color=OI["vermilion"],
            markeredgecolor="k", markeredgewidth=0.4, zorder=7)
    # detector block on the bus
    ax.add_patch(Rectangle((10.3, 0.65), 0.55, 0.42, fc=DET, ec="k",
                           lw=0.5))
    # labels
    ax.text(1.7, 1.13, "diamond, 200 nm", fontsize=7, ha="center",
            va="center")
    ax.plot([1.7, 1.7], [0.95, 0.665], color="#999999", lw=0.7)
    ax.text(8.7, 1.13, "TFLN film, 190 nm", fontsize=7, ha="center",
            va="center")
    ax.plot([8.8, 8.9], [0.95, 0.655], color="#999999", lw=0.7)
    ax.text(0.22, 0.23, "SiO$_2$", fontsize=7, color="#666666",
            va="center")
    ax.text(10.15, 0.23, "side view, not to scale", fontsize=6,
            color="#999999", ha="right", va="center")


def panel_zoom(ax):
    """Magnified plan view of the cavity region."""
    ax.set_xlim(-0.4, 13.4); ax.set_ylim(-0.35, 2.65)
    ax.set_aspect("equal"); ax.axis("off")
    A = 0.55                      # mirror period (scene units)
    yc, hw = 1.05, 0.40           # beam centre, half width (w = 1.44a)
    ax.add_patch(Rectangle((-0.4, -0.35), 13.8, 3.0, fc=SUB, ec="none"))
    ax.add_patch(Rectangle((0.1, yc-hw), 12.4, 2*hw, fc=DIA, ec="k",
                           lw=0.5, zorder=4))
    # gaps: 3 mirror | quadratic taper down p1..p8 | up p7..p1 | 2 output
    pd = [A*(1 - 0.14*(i/8.0)**2) for i in range(1, 9)]
    gaps = [A]*3 + pd + pd[-2::-1] + [A]*2
    xs = [0.62]
    for g in gaps:
        xs.append(xs[-1] + g)
    for xh, g in zip(xs, [A] + gaps):
        ax.add_patch(plt.Circle((xh, yc), 0.30*g, fc="white", ec="k",
                                lw=0.35, zorder=5))
    x_cav = 0.62 + 3*A + sum(pd) - pd[-1]/2
    # mode glow + emitter
    for rr, aa in ((1.0, 0.30), (0.55, 0.5)):
        ax.add_patch(Ellipse((x_cav, yc), 2.1*rr, 0.62*rr,
                             fc=OI["orange"], ec="none", alpha=aa,
                             zorder=6))
    ax.plot([x_cav], [yc], marker="*", ms=9, color=OI["vermilion"],
            markeredgecolor="k", markeredgewidth=0.4, zorder=7)
    # continuation dots at the cuts
    ax.text(-0.05, yc, "$\\cdots$", fontsize=9, ha="right", va="center")
    ax.text(11.82, yc, "$\\cdots$", fontsize=9, ha="left", va="center")
    # emission toward the taper and bus
    ax.annotate("", xy=(13.35, yc), xytext=(12.58, yc),
                arrowprops=dict(arrowstyle="-|>", lw=1.2,
                                color=OI["green"]), zorder=7,
                annotation_clip=False)
    ax.text(13.0, yc + 0.45, "to bus", fontsize=6.8, color=OI["green"],
            ha="right")
    # region brackets above the beam
    def bracket(x0, x1, label, ty=2.30):
        yb = 1.80
        ax.plot([x0, x1], [yb, yb], color="k", lw=0.7)
        ax.plot([x0, x0], [yb-0.12, yb], color="k", lw=0.7)
        ax.plot([x1, x1], [yb-0.12, yb], color="k", lw=0.7)
        ax.text((x0+x1)/2, ty, label, fontsize=6.8, ha="center",
                va="top")
    bracket(0.1, 2.30, "mirror, $N = 10$")
    bracket(2.42, 9.86, "quadratic taper, $a \\rightarrow 0.86a$")
    bracket(9.98, 11.35, "output, $Q_L \\approx 500$")
    # lattice-constant dimension between two mirror holes
    for xe in (xs[1], xs[2]):
        ax.plot([xe, xe], [yc - 0.16, 0.38], color="#777777", lw=0.5)
    ax.plot([xs[1], xs[2]], [0.44, 0.44], color="#333333", lw=0.8)
    ax.text((xs[1]+xs[2])/2, 0.26, "$a = 194$ nm", fontsize=6.8,
            ha="center", va="top")
    # hole radius
    _lead(ax, (4.1, 0.02), (xs[3] + 0.12, yc - 0.16), "$r = 0.30a$",
          fs=7, va="top")
    # emitter label
    _lead(ax, (7.9, 0.02), (x_cav + 0.12, yc - 0.24), "SnV$^-$ spin",
          fs=7, color=OI["vermilion"], va="top")


def panel_b(ax):
    ax.axis("off")
    x0, x1 = 0.14, 0.64
    yC = 0.88; y1 = 0.13; y2 = 0.29
    # levels: the |up> line starts to the right so the C1 arrow does not
    # cross it
    ax.plot([x0, x1], [yC, yC], color="k", lw=1.1)
    ax.plot([x0, x1], [y1, y1], color="k", lw=1.1)
    ax.plot([0.345, x1], [y2, y2], color="k", lw=1.1)
    for y, lab in ((yC, "$|C\\rangle$"), (y1, "$|\\!\\downarrow\\rangle$"),
                   (y2, "$|\\!\\uparrow\\rangle$")):
        ax.text(x1 + 0.02, y, lab, fontsize=8, va="center")
    # qubit splitting bracket next to the |up> line start
    ax.annotate("", xy=(0.44, y2), xytext=(0.44, y1),
                arrowprops=dict(arrowstyle="<->", lw=0.7,
                                mutation_scale=6))
    ax.text(0.46, (y1+y2)/2, "$\\omega_q$", fontsize=7, ha="left",
            va="center")
    # spin-conserving C1 (cavity-enhanced), clear of the |up> line
    ax.add_patch(FancyArrowPatch((0.195, y1), (0.195, yC),
                 arrowstyle="<|-|>", mutation_scale=8, lw=1.4,
                 color=OI["vermilion"]))
    ax.text(0.175, 0.56, "C$_1$", fontsize=7.5, color=OI["vermilion"],
            ha="right")
    # cavity line profile beside C1
    yy = np.linspace(0.53, 0.84, 120)
    lor = 0.075 / (1 + ((yy - 0.67) / 0.085) ** 2)
    ax.plot(0.395 + lor, yy, color=OI["green"], lw=1.2)
    ax.text(0.36, 0.41, "cavity", fontsize=6.5, ha="center",
            color=OI["green"])
    ax.text(0.36, 0.32, "$F_C = 19$", fontsize=6.5, ha="center",
            color=OI["green"])
    # weak spin-flipping C2
    ax.add_patch(FancyArrowPatch((0.54, y2), (0.54, yC), arrowstyle="-|>",
                 mutation_scale=7, lw=1.0, color=OI["grey"],
                 linestyle=(0, (3, 2))))
    ax.text(0.565, 0.50, "C$_2$ (weak,\n$\\propto 1/\\Lambda_0$)",
            fontsize=6.8, color=OI["grey"])
    ax.text(0.5, 0.0,
            "$\\gamma_{\\rm cav}=\\zeta\\gamma \\ll \\omega_q$:"
            " readout stays spin-selective",
            fontsize=6.8, ha="center")
    ax.set_xlim(0, 1); ax.set_ylim(-0.03, 1)


def panel_c(ax):
    L0 = np.array(frl["Lambda0"])
    ax.semilogx(L0, 100*np.array(frl["cavity"]), color=OI["green"], lw=1.5)
    ax.fill_between(L0, 100*np.array(frl["confocal_002"]),
                    100*np.array(frl["confocal_004"]), color=OI["grey"],
                    alpha=0.45, lw=0)
    ax.plot([2244], [87.4], "o", ms=4, color=OI["vermilion"])
    ax.text(60, 90.0, "this design", color=OI["green"], fontsize=6.5)
    ax.text(700, 63.5, "confocal", color="#777777", fontsize=6.5,
            ha="center")
    ax.text(2244, 80.5, "measured", color=OI["vermilion"], fontsize=6.5,
            ha="center", va="top")
    ax.set_xlabel("bare cyclicity $\\Lambda_0$", fontsize=8)
    ax.set_ylabel("$F_r$ (%)", fontsize=8)
    ax.set_ylim(48, 103)
    ax.annotate("98.5% in 0.1 µs", xy=(2244, 98.5), xytext=(9, 68),
                fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.6))



fig = plt.figure(figsize=(DOUBLE, 2.40))
gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1], wspace=0.20,
                      left=0.005, right=0.985, top=0.945, bottom=0.115)
gsl = gs[0, 0].subgridspec(3, 1, height_ratios=[1.0, 0.44, 0.92],
                           hspace=0.10)
gsr = gs[0, 1].subgridspec(2, 1, height_ratios=[1.05, 1], hspace=0.42)
axa1 = fig.add_subplot(gsl[0]); panel_top(axa1); axa1.set_anchor("W")
axa2 = fig.add_subplot(gsl[1]); panel_side(axa2); axa2.set_anchor("W")
axa3 = fig.add_subplot(gsl[2]); panel_zoom(axa3); axa3.set_anchor("W")
axa1.text(-0.005, 1.06, "(a)", transform=axa1.transAxes,
          fontweight="bold", fontsize=9, va="bottom")
axb = fig.add_subplot(gsr[0]); panel_b(axb); panel_label(axb, "(b)", dx=-0.06)
axc = fig.add_subplot(gsr[1]); panel_c(axc); panel_label(axc, "(c)", dx=-0.28)
fig.savefig(HERE / "fig1_device.pdf")
fig.savefig(HERE / "fig1_device.png")
print("saved fig1_device")
