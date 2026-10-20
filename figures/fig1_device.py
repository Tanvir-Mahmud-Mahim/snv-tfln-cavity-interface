"""Fig. 1: device concept (2D top + side views with a magnified cavity
region), SnV- level structure with all readout transitions, and the
optical spectrum around the driven A1 line at the reference operating
point (revision 1)."""
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
    # gaps: 3 mirror | quadratic taper down p1..p8 | up p7..p1 | N_out
    # output mirror periods of the design (FDTD design point)
    n_out = json.load(open(RES / "cavity_loaded.json"))["n_out"]
    pd = [A*(1 - 0.14*(i/8.0)**2) for i in range(1, 9)]
    gaps = [A]*3 + pd + pd[-2::-1] + [A]*n_out
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
    def bracket(x0, x1, label, ty=2.30, ha="center"):
        yb = 1.80
        ax.plot([x0, x1], [yb, yb], color="k", lw=0.7)
        ax.plot([x0, x0], [yb-0.12, yb], color="k", lw=0.7)
        ax.plot([x1, x1], [yb-0.12, yb], color="k", lw=0.7)
        ax.text((x0+x1)/2 if ha == "center" else x0, ty, label,
                fontsize=6.8, ha=ha, va="top")
    bracket(0.1, 2.30, "mirror, $N = 10$")
    x_t1 = xs[-1 - n_out] + 0.30          # end of the output-side taper
    bracket(2.42, x_t1, "quadratic taper, $a \\rightarrow 0.86a$")
    QL = json.load(open(RES / "cavity_loaded.json"))["Q_L"]
    bracket(x_t1 + 0.12, 11.6, f"output, $Q_L \\approx {QL:.0f}$", ha="left")
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
    """Level structure of the SnV- with all transitions used by the
    readout model (not to scale; energy axis broken twice)."""
    ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    xl, xr = 0.13, 0.72
    yD, yU, yUG = 0.06, 0.17, 0.33          # |dn>, |up>, upper ground
    yA, yB, yUE = 0.66, 0.77, 0.95          # |A>, |B>, upper excited
    for y in (yD, yU, yA, yB):
        ax.plot([xl, xr], [y, y], color="k", lw=1.0)
    for y in (yUG, yUE):
        ax.plot([xl, xr], [y, y], color="#777777", lw=0.8, ls=(0, (4, 2)))
    lab = {yD: "$|\\!\\downarrow\\rangle$", yU: "$|\\!\\uparrow\\rangle$",
           yA: "$|A\\rangle$", yB: "$|B\\rangle$"}
    for y, t in lab.items():
        ax.text(xr + 0.02, y, t, fontsize=7.5, va="center")
    ax.text(xr + 0.02, yUG, "upper\nground", fontsize=6.5, va="center",
            color="#555555", linespacing=0.95)
    ax.text(xr + 0.02, yUE, "upper\nexcited", fontsize=6.5, va="center",
            color="#555555", linespacing=0.95)
    # splittings (left side): omega_q and Delta_g stacked, omega_e
    for (ya, yb, t, c) in ((yD, yU, "$\\omega_q$", "k"),
                           (yA, yB, "$\\omega_e$", "k"),
                           (yU, yUG, "$\\Delta_g$", "#555555")):
        ax.annotate("", xy=(0.075, yb), xytext=(0.075, ya),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color=c,
                                    mutation_scale=5, shrinkA=0, shrinkB=0))
        ax.text(0.06, (ya + yb) / 2, t, fontsize=7, ha="right", va="center",
                color=c)
    # transitions
    def arrow(x, y0, y1, color, lw, ls="-", style="-|>"):
        ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle=style,
                     mutation_scale=6, lw=lw, color=color, linestyle=ls,
                     shrinkA=0, shrinkB=0))
    arrow(0.20, yD, yA, OI["vermilion"], 1.5, style="<|-|>")
    arrow(0.33, yU, yB, OI["blue"], 1.1)
    arrow(0.45, yU, yA, OI["grey"], 0.8, ls=(0, (2, 1.5)))
    arrow(0.55, yD, yB, OI["grey"], 0.8, ls=(0, (2, 1.5)))
    arrow(0.65, yA, yUG, "#777777", 0.7, ls=(0, (1, 1.2)))
    for x, t, c in ((0.20, "A1", OI["vermilion"]), (0.33, "B2", OI["blue"]),
                    (0.45, "A2", "#666666"), (0.55, "B1", "#666666"),
                    (0.65, "D", "#666666")):
        ax.text(x + 0.012, 0.50 if t != "D" else 0.55, t, fontsize=6.8,
                color=c, ha="left", va="center")


def panel_c(ax):
    """Optical spectrum around A1 at the reference operating point."""
    rr = json.load(open(RES / "readout_rev.json"))["reference"]
    lv = rr["levels"]
    g_cav = rr["budget"]["gamma_GHz"]
    g0 = 1.0 / (2 * np.pi * 4.5)
    S = np.array(lv["S_rel"])
    pos = {"A1": 0.0, "B2": lv["dA1B2"], "A2": lv["dA1A2"],
           "B1": lv["dA1B1"]}
    stren = {"A1": S[0][0], "B2": S[1][1], "A2": S[1][0], "B1": S[0][1]}
    col = {"A1": OI["vermilion"], "B2": OI["blue"], "A2": OI["grey"],
           "B1": OI["grey"]}
    x = np.linspace(-5.2, 5.6, 4000)
    for k in ("A1", "B2", "A2", "B1"):
        L = stren[k] / (1 + (2 * (x - pos[k]) / g_cav) ** 2)
        ax.semilogy(x, L, color=col[k], lw=1.1 if k in ("A1", "B2") else 0.8)
    Lb = stren["B2"] / (1 + (2 * (x - pos["B2"]) / g0) ** 2)
    ax.semilogy(x, Lb, color=OI["blue"], lw=0.8, ls=(0, (1.5, 1.2)))
    ax.axvline(0.0, color="k", lw=0.5, ls=(0, (3, 2)))
    leak = stren["B2"] / (1 + (2 * pos["B2"] / g_cav) ** 2)
    ax.plot([0], [leak], "o", ms=3, color=OI["blue"])
    ax.set_ylim(2e-6, 4)
    ax.set_xlim(-5.2, 5.6)
    ax.set_xlabel("detuning from A1 (GHz)")
    ax.set_ylabel("rel. excitation")
    ax.set_yticks([1e-4, 1e-2, 1])
    ax.set_xticks([-4, -2, 0, 2, 4])
    for k in ("A2", "B1"):
        ax.text(pos[k], 4e-3, k, fontsize=6.5, ha="center", va="bottom",
                color="#666666")
    ax.text(1.15, 1.2, "B2", fontsize=6.5, ha="left", va="center",
            color=OI["blue"])
    ax.text(-0.45, 1.2, "A1", fontsize=6.5, ha="right", va="center",
            color=OI["vermilion"])
    ax.annotate("dark-state\nleak", xy=(-0.05, leak), xytext=(-3.3, 0.25),
                fontsize=6.3, ha="center", va="center", linespacing=0.95,
                arrowprops=dict(arrowstyle="->", lw=0.5, shrinkB=2))


fig = plt.figure(figsize=(DOUBLE, 2.40))
gs = fig.add_gridspec(1, 2, width_ratios=[1.8, 1], wspace=0.20,
                      left=0.005, right=0.985, top=0.945, bottom=0.115)
gsl = gs[0, 0].subgridspec(3, 1, height_ratios=[1.0, 0.44, 0.92],
                           hspace=0.10)
gsr = gs[0, 1].subgridspec(2, 1, height_ratios=[1.45, 0.75], hspace=0.38)
axa1 = fig.add_subplot(gsl[0]); panel_top(axa1); axa1.set_anchor("W")
axa2 = fig.add_subplot(gsl[1]); panel_side(axa2); axa2.set_anchor("W")
axa3 = fig.add_subplot(gsl[2]); panel_zoom(axa3); axa3.set_anchor("W")
axa1.text(-0.005, 1.06, "(a)", transform=axa1.transAxes,
          fontweight="bold", fontsize=9, va="bottom")
axb = fig.add_subplot(gsr[0]); panel_b(axb); panel_label(axb, "(b)", dx=-0.02, dy=1.03)
axc = fig.add_subplot(gsr[1]); panel_c(axc)
fig.canvas.draw()
_p = axc.get_position()
fig.text(_p.x0 - 0.085, _p.y1 + 0.035, "(c)", fontweight="bold", fontsize=9,
         ha="left", va="bottom")
fig.savefig(HERE / "fig1_device.pdf")
fig.savefig(HERE / "fig1_device.png")
print("saved fig1_device")
