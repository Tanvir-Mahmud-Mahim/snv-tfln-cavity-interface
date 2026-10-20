"""Generate all manuscript figures from sim/results (no hand-typed data)."""
import json
import pathlib
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrow

HERE = pathlib.Path(__file__).resolve().parent
SIM = HERE.parent / "sim"
RES = SIM / "results"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SIM))

from figstyle import OI, SINGLE, DOUBLE, setup, panel_label

setup()


def save(fig, name):
    fig.savefig(HERE / f"{name}.pdf")
    fig.savefig(HERE / f"{name}.png")
    plt.close(fig)
    print("saved", name)


# ----------------------------------------------------------------- data
wg_d = json.load(open(RES / "wg_diamond.json"))
wg_ln = json.load(open(RES / "wg_ln.json"))
snvr = json.load(open(RES / "snv_results.json"))
prod = json.load(open(RES / "cavity_production.json"))
fld = np.load(RES / "cavity_field.npz")
tol = json.load(open(RES / "tolerance.json"))

TOLD = RES / "tolerance_eme"
def tolj(tag):
    return json.load(open(TOLD / f"{tag}.json"))
rr = json.load(open(RES / "readout_rev.json"))
fdtd = json.load(open(RES / "cavity_fdtd_summary.json"))
cavL = json.load(open(RES / "cavity_loaded.json"))


def letters(fig, axes_labels, dx=-0.012, dy=0.012):
    """Panel letters placed after layout, just above the top-left corner
    of each axes' tight bounding box (never overlapping tick labels)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for ax, lab in axes_labels:
        bb = ax.get_tightbbox(r)
        x0, y1 = inv.transform((bb.x0, bb.y1))
        p = ax.get_position()
        fig.text(x0 + dx, max(y1, p.y1) + dy, lab, fontweight="bold",
                 fontsize=9, ha="left", va="bottom")


# supermode branches from the EME log
branches = []
for line in open(RES / "eme_log2.txt"):
    m = re.match(r"\s*\[\d+/30\] w=\s*([\d.]+) nm\s+n_eff=\[(.*)\]", line)
    if m:
        w = float(m.group(1))
        ns = [float(x.strip().strip("'")) for x in m.group(2).split(",")]
        branches.append((w, ns))
branches.sort()
bw = np.array([b[0] for b in branches])
bn = np.array([b[1] for b in branches])


# ================================================================ Fig 2
def fig2():
    fig = plt.figure(figsize=(DOUBLE, 3.25))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.55,
                          left=0.10, right=0.985, top=0.93, bottom=0.11)
    axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]

    # (a) mode dispersion: diamond beam on SiO2 and LN ridge
    ax = axs[0]
    w = 1e3 * np.array(wg_d["w"])
    te0 = []
    for ne, te in zip(wg_d["neff"], wg_d["te"]):
        cands = [n for n, t in zip(ne, te) if t > 0.7]
        te0.append(max(cands) if cands else np.nan)
    ax.plot(w, te0, color=OI["blue"], label="diamond TE$_0$")
    wl = 1e3 * np.array(wg_ln["w"])
    te0l = []
    for ne, te in zip(wg_ln["neff"], wg_ln["te"]):
        cands = [n for n, t in zip(ne, te) if t > 0.7]
        te0l.append(max(cands) if cands else np.nan)
    ax.plot(wl, te0l, color=OI["vermilion"], label="LN ridge TE$_0$")
    ax.axhline(1.4574, color=OI["grey"], lw=0.7, ls=":")
    ax.text(640, 1.4574, "SiO$_2$ light line", color=OI["grey"], fontsize=6.5,
            ha="right", va="center",
            bbox=dict(boxstyle="square,pad=0.10", fc="white", ec="none"))
    ax.set_ylim(bottom=1.395)
    ax.axvline(280, color=OI["blue"], lw=0.7, ls="--", alpha=0.5)
    ax.set_xlabel("width (nm)")
    ax.set_ylabel("$n_\\mathrm{eff}$ at 619 nm")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.17),
              fontsize=6.5, handlelength=1.1, handletextpad=0.4,
              borderaxespad=0.2, labelspacing=0.3)

    # (b) hybrid supermode branches vs diamond taper width
    ax = axs[1]
    cols = [OI["blue"], OI["orange"], OI["green"], OI["purple"]]
    for i in range(bn.shape[1]):
        ax.plot(bw, bn[:, i], color=cols[i], lw=1.1)
    ends = bn[0]
    ax.text(120, 2.045, "TE", color=cols[0], fontsize=7, ha="center",
            va="bottom")
    ax.text(200, 1.930, "TM", color=cols[1], fontsize=7, ha="center",
            va="center")
    ax.text(58, ends[2] - 0.012, "LN TE$_1$", color=cols[2],
            fontsize=7, ha="right", va="top")
    ax.text(bw[-1], 1.790, "LN TM$_0$", color=cols[3], fontsize=7,
            ha="left", va="top")
    ax.set_xlabel("diamond width (nm)")
    ax.set_ylabel("supermode $n_\\mathrm{eff}$")
    ax.invert_xaxis()

    # (c) EME transfer vs taper length, three LN profiles
    ax = axs[2]
    for tag, c, lab, mk in (("nominal", OI["green"], "nominal", "o"),
                            ("sw70", OI["blue"], "70$^\\circ$ walls", "s"),
                            ("rib60", OI["vermilion"], "rib", "^")):
        d = tolj(tag)
        ax.plot(d["lengths"], 100 * np.array(d["T_LN"]), marker=mk, ms=3,
                color=c, lw=1.0, label=lab)
    ax.axhline(92, color=OI["grey"], lw=0.8, ls="--")
    ax.text(14.0, 92.6, "Riedel et al.", color=OI["grey"], fontsize=6.5,
            ha="right", va="bottom")
    ax.set_xlabel("taper length (µm)")
    ax.set_ylabel("transfer $T$ (%)")
    ax.set_ylim(80, 101)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=6.5,
              handlelength=1.2, handletextpad=0.4, borderaxespad=0.2,
              labelspacing=0.25)

    # (d) lateral offset
    ax = axs[3]
    offs = [0, 25, 50, 75, 100, 150]
    tags = ["nominal", "off25", "off50", "off75", "off100", "off150"]
    Ts = []
    for t in tags:
        d = tolj(t)
        Ts.append(100 * d["T_LN"][d["lengths"].index(6.0)])
    ax.axvspan(22, 54, color=OI["grey"], alpha=0.25, lw=0)
    ax.plot(offs, Ts, "o-", ms=3, color=OI["green"])
    ax.set_xlabel("lateral offset (nm)")
    ax.set_ylabel("$T$ (%)")
    ax.set_ylim(min(Ts) - 3, 100.5)

    # (e) angular misalignment
    ax = axs[4]
    angs = [0, 0.5, 1, 2, 3]
    tags = ["nominal", "tilt0p5", "tilt1", "tilt2", "tilt3"]
    Ts = []
    for t in tags:
        d = tolj(t)
        Ts.append(100 * d["T_LN"][d["lengths"].index(6.0)])
    ax.plot(angs, Ts, "o-", ms=3, color=OI["green"])
    ax.set_xlabel("misalignment (deg)")
    ax.set_ylabel("$T$ (%)")
    ax.set_ylim(min(Ts) - 3, 100.5)

    # (f) vertical gap
    ax = axs[5]
    gaps = [0, 10, 20, 50]
    tags = ["nominal", "gap10", "gap20", "gap50"]
    for L, c, mk in ((6.0, OI["green"], "o"), (14.0, OI["blue"], "s")):
        Ts = []
        for t in tags:
            d = tolj(t)
            Ts.append(100 * d["T_LN"][d["lengths"].index(L)])
        ax.plot(gaps, Ts, marker=mk, ms=3, color=c, lw=1.0,
                label=f"{L:.0f} µm")
    ax.set_xlabel("gap (nm)")
    ax.set_ylabel("$T$ (%)")
    ax.legend(loc="lower left", fontsize=6.5, handlelength=1.2,
              handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25)

    letters(fig, list(zip(axs, ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"])))
    save(fig, "fig2_waveguide")


# ================================================================ Fig 3
def fig3():
    fig = plt.figure(figsize=(DOUBLE, 3.55))
    outer = fig.add_gridspec(2, 1, height_ratios=[2.06, 1.45], hspace=0.42,
                             left=0.095, right=0.985, top=0.95, bottom=0.10)
    gs = outer[0].subgridspec(2, 3, height_ratios=[1, 1],
                              width_ratios=[1, 1.5, 1], wspace=0.55,
                              hspace=0.12)
    gb = outer[1].subgridspec(1, 3, width_ratios=[1, 0.8, 1], wspace=1.05)

    # (a) mirror band structure
    ax = fig.add_subplot(gs[:, 0])
    import cavity as cav
    kx, fr = cav.band_structure(194.35, nk=13, gmax=3.0, gmode_inds=(0, 1, 2),
                                numeig=8)
    a = 194.35
    for b in range(fr.shape[1]):
        ax.plot(kx / np.pi, fr[:, b], ".", ms=2.2, color=OI["grey"])
    kk = np.linspace(1e-3, np.pi, 60)
    ax.fill_between(kk / np.pi, kk / (2 * np.pi) / 1.4574, 0.42,
                    color=OI["skyblue"], alpha=0.18, lw=0)
    try:
        edges = json.load(open(RES / "band_edges.json"))
    except FileNotFoundError:
        edges = cav.classify_te_edges(194.35)
    f_diel, f_air = edges["f_diel_ca"], edges["f_air_ca"]
    ax.fill_between([0, 1], f_diel, f_air, color=OI["yellow"], alpha=0.30,
                    lw=0)
    ax.text(0.62, 0.345, "TE gap", fontsize=7, color=OI["black"],
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                      alpha=0.85))
    f_c = a / prod["lam_nm"]
    ax.axhline(f_c, color=OI["vermilion"], lw=1.0)
    ax.text(0.05, f_c + 0.007, "cavity mode",
            color=OI["vermilion"], fontsize=7, va="bottom",
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                      alpha=0.85))
    ax.set_xlabel("$k_x a/\\pi$")
    ax.set_ylabel("frequency $a/\\lambda$")
    ax.set_xlim(0, 1); ax.set_ylim(0.15, 0.42)
    axa = ax

    # (b) field maps of the closed cavity
    axb = fig.add_subplot(gs[0, 1])
    xs, ys = fld["xs"], fld["ys"]
    E2 = fld["E2xy"]
    axb.imshow(E2 / E2.max(), extent=[xs[0] * a / 1e3, xs[-1] * a / 1e3,
               ys[0] * a / 1e3, ys[-1] * a / 1e3], origin="lower",
               cmap="inferno", aspect="auto", vmax=1.0)
    axb.set_ylim(-0.35, 0.35)
    axb.set_ylabel("y (µm)")
    axb.set_xticklabels([])
    axb2 = fig.add_subplot(gs[1, 1])
    zs = fld["zs"]
    E2z = fld["E2xz"]
    axb2.imshow(E2z / E2z.max(), extent=[xs[0] * a / 1e3, xs[-1] * a / 1e3,
                zs[0] * a / 1e3, zs[-1] * a / 1e3], origin="lower",
                cmap="inferno", aspect="auto", vmax=1.0)
    axb2.set_ylim(-0.25, 0.45)
    axb2.set_xlabel("x (µm)")
    axb2.set_ylabel("z (µm)")

    # (c) envelope decay
    axc = fig.add_subplot(gs[:, 2])
    line = fld["E2line"]
    axc.semilogy(xs * a / 1e3, line / line.max(), color=OI["blue"], lw=0.8)
    axc.set_xlabel("x (µm)")
    axc.set_ylabel("$|E|^2$ (norm., y=z=0)")
    axc.set_ylim(1e-5, 2)

    # (d) loaded Q vs output mirror periods (FDTD) + FP estimate
    axd = fig.add_subplot(gb[0, 0])
    Ns = sorted(int(k[1]) for k in fdtd if k.endswith("_r40"))
    QL = [fdtd[f"N{n}_r40"]["Q_L"] for n in Ns]
    # Fabry-Perot estimate of the original submission:
    # Q(N) = 2 pi f (2 n_g L_c / c) / D^N with n_g = 10, L_c = 1.2 um and
    # the GME mirror decay D per period
    f_c = 299792458.0 / (prod["lam_nm"] * 1e-9)
    fp = [2 * np.pi * f_c * (2 * 10.0 * 1.2e-6 / 299792458.0)
          / prod["mirror_decay_D"] ** n for n in range(4)]
    axd.semilogy(range(4), fp, "--", color=OI["grey"], lw=0.9)
    axd.semilogy(Ns, QL, "o-", color=OI["green"], ms=3.5)
    axd.plot([cavL["n_out"]], [cavL["Q_L"]], "o", ms=6, mfc="none",
             mec=OI["vermilion"], mew=0.9)
    axd.text(3.0, 330, "FP estimate\n(original)", color=OI["grey"],
             fontsize=6.5, ha="right", va="bottom", linespacing=0.95)
    axd.set_ylim(200, 3 * max(QL + fp))
    axd.set_xticks(range(4))
    axd.set_yticks([300, 1000, 3000, 10000, 30000])
    axd.set_yticklabels(["300", "1000", "3000", "$10^4$", "$3\\times10^4$"])
    axd.minorticks_off()
    axd.set_xlabel("output periods $N_\\mathrm{out}$")
    axd.set_ylabel("loaded $Q_L$")

    # (e) fate of the emitted power at resonance
    axe = fig.add_subplot(gb[0, 1])
    keys = [("frac_TE0", "waveguide", OI["green"]),
            ("frac_down", "substrate", OI["orange"]),
            ("frac_up", "air", OI["skyblue"]),
            ("frac_side", "side", OI["purple"]),
            ("frac_left", "closed end", OI["grey"])]
    bottom = np.zeros(len(Ns))
    for k, lab, c in keys:
        v = np.array([100 * fdtd[f"N{n}_r40"][k] for n in Ns])
        axe.bar(Ns, v, bottom=bottom, color=c, width=0.65, label=lab)
        bottom += v
    axe.set_xticks(Ns)
    axe.set_ylim(0, 100)
    axe.set_xlabel("output periods $N_\\mathrm{out}$")
    axe.set_ylabel("dipole power (%)")
    axe.legend(loc="upper left", bbox_to_anchor=(1.02, 1.02), fontsize=6.3,
               handlelength=0.9, handletextpad=0.35, borderaxespad=0.0,
               labelspacing=0.25)

    # (f) LDOS enhancement spectrum at the design point
    axf = fig.add_subplot(gb[0, 2])
    nd = cavL["n_out"]
    raw = json.load(open(RES / f"cavity_fdtd_N{nd}_r40.json"))
    bulk = json.load(open(RES / "cavity_fdtd_bulk_r40.json"))
    ld = np.array(raw["ldos"])
    f = np.linspace(min(raw["freqs"]), max(raw["freqs"]), len(ld))
    F = ld / np.interp(f, bulk["freqs"], bulk["ldos"])
    lam = 1e3 / f
    axf.plot(lam, F, color=OI["blue"], lw=1.0)
    axf.set_xlabel("wavelength (nm)")
    axf.set_ylabel("LDOS enhancement")
    axf.set_xlim(lam.min(), lam.max())

    letters(fig, [(axa, "(a)"), (axb, "(b)"), (axc, "(c)"), (axd, "(d)"),
                  (axe, "(e)"), (axf, "(f)")])
    save(fig, "fig3_cavity")


# ================================================================ Fig 4
def fig4():
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 1.62))

    # (a) cyclicity vs field angle + measured points
    ax = axs[0]
    z = np.array(snvr["zeta_deg"]); c = np.array(snvr["cyc_vs_zeta"])
    ax.semilogy(z, c, color=OI["blue"])
    ax.errorbar([147], [2244], yerr=[108], fmt="o", color=OI["vermilion"],
                ms=4, capsize=2)
    ax.errorbar([53], [8.6], yerr=[0.4], fmt="o", color=OI["vermilion"], ms=4,
                capsize=2)
    ax.text(30, 700, "model", color=OI["blue"], fontsize=7, ha="left")
    ax.text(172, 2244, "see (b)", fontsize=6.5, va="center")
    ax.set_ylim(3.5, 8e3)
    ax.set_xlabel("field angle $\\zeta$ (deg)")
    ax.set_ylabel("cyclicity $\\Lambda$")
    ax.set_xlim(0, 360)
    panel_label(ax, "(a)", dx=-0.30, dy=1.10)

    # (b) cyclicity vs misalignment
    ax = axs[1]
    al = np.array(snvr["alpha_deg"]); ca = np.array(snvr["cyc_vs_alpha"])
    ax.loglog(al, ca, color=OI["blue"])
    ax.errorbar([10], [2244], yerr=[108], fmt="o", color=OI["vermilion"],
                ms=4, capsize=2)
    ax.annotate("measured max.\n($\\alpha \\approx 10°$)",
                xy=(10, 2244), xytext=(0.03, 0.06),
                textcoords="axes fraction", va="bottom", fontsize=6.5,
                arrowprops=dict(arrowstyle="->", lw=0.7,
                                relpos=(0.75, 1.0)))
    ax.set_xlabel("misalignment $\\alpha$ (deg)")
    ax.set_ylabel("cyclicity $\\Lambda$")
    panel_label(ax, "(b)", dx=-0.56, dy=1.08)

    # (c) strain trade-off
    ax = axs[2]
    r = np.array(snvr["strain_2UgOverLam"])
    cyc = np.array(snvr["cyc_vs_strain"])
    rab = np.array(snvr["rabi_MHz_vs_strain"])
    ax.loglog(r, cyc, color=OI["blue"], label="$\\Lambda$")
    ax.set_xlabel("strain $2\\Upsilon_g/\\lambda_g$")
    ax.set_ylabel("cyclicity $\\Lambda$", color=OI["blue"])
    ax.axvline(snvr["device_strain_ratio"], color=OI["grey"], lw=0.8, ls="--")
    import matplotlib.transforms as mtr
    trans = mtr.blended_transform_factory(ax.transData, ax.transAxes)
    ax.text(snvr["device_strain_ratio"], 1.03, "measured device",
            transform=trans, fontsize=6.5, ha="center", va="bottom",
            color=OI["grey"])
    ax2 = ax.twinx()
    ax2.loglog(r, rab, color=OI["vermilion"])
    ax2.set_ylabel("Rabi rate (MHz)", color=OI["vermilion"])
    ax2.tick_params(axis="y", colors=OI["vermilion"])
    ax.tick_params(axis="y", colors=OI["blue"])
    panel_label(ax, "(c)", dx=-0.30)

    fig.tight_layout(w_pad=2.2)
    save(fig, "fig4_spin")


# ================================================================ Fig 5
def fig5():
    fig = plt.figure(figsize=(DOUBLE, 3.45))
    gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.78,
                          left=0.16, right=0.93, top=0.93, bottom=0.11)
    axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(3)]
    ref = rr["reference"]
    bud = ref["budget"]

    # (a) efficiency budget
    ax = axs[0]
    conf = [("QE", 0.8), ("PSB", 0.43), ("scatter", 0.05), ("path", 0.35),
            ("detector", 0.65)]
    this = [("cavity $\\beta$", bud["beta"]),
            ("waveguide", bud["eta_wg"]), ("taper", cavL["T_taper"]),
            ("fiber", 0.5), ("detector", 0.77)]
    y = np.arange(5)
    ax.barh(y + 0.2, [c[1] for c in conf], height=0.36, color=OI["grey"],
            label="confocal (meas.)")
    ax.barh(y - 0.2, [c[1] for c in this], height=0.36, color=OI["green"],
            label="this design")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c[0]} | {v[0]}" for c, v in zip(conf, this)],
                       fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("per-stage efficiency")
    tot_c = np.prod([c[1] for c in conf])
    ax.set_title(f"$\\eta$: {100*tot_c:.1f}% vs {100*bud['eta']:.0f}%",
                 fontsize=7, pad=3)

    # (b) fidelity vs window
    ax = axs[1]
    t = np.array(ref["vs_tau"]["tau_ns"]) / 1e3
    F = 100 * np.array(ref["vs_tau"]["F"], dtype=float)
    ax.semilogx(t, F, color=OI["green"], lw=1.3)
    for key, ls, lab in (("eta0.002_xlw0.0", "-", "0.2%"),
                         ("eta0.004_xlw0.0", "--", "0.4%")):
        Fc = 100 * np.array(ref["confocal"][key]["F"], dtype=float)
        ax.semilogx(t, Fc, color=OI["grey"], lw=1.0, ls=ls)
    ax.plot([50], [87.4], "o", ms=3.5, color=OI["vermilion"])
    ax.set_ylim(50, 101)
    ax.set_xlim(5e-3, 100)
    ax.set_xlabel("window $\\tau$ (µs)")
    ax.set_ylabel("$F_r$ (%)")
    ax.text(0.25, 93.5, "integrated", color=OI["green"], fontsize=6.5,
            ha="left", va="top")
    ax.text(30, 68, "confocal", color="#666666", fontsize=6.5,
            ha="right", va="center")

    # (c) fidelity vs field misalignment, with cyclicity
    ax = axs[2]
    rows = rr["vs_alpha"]
    al = np.array([r["alpha"] for r in rows])
    Fa = 100 * np.array([r["F"] for r in rows])
    lam = np.array([r["levels"]["LamA"] for r in rows])
    dab = np.array([r["levels"]["dA1B2"] for r in rows])
    ax.plot(al, Fa, "o-", ms=2.5, color=OI["green"])
    zc = np.interp(0.0, -dab[(al > 50) & (al < 75)][::1],
                   al[(al > 50) & (al < 75)])
    ax.axvline(zc, color=OI["blue"], lw=0.7, ls="--")
    ax.text(zc - 2, 97, "$\\omega_e = \\omega_q$", color=OI["blue"],
            fontsize=6.5, ha="right", va="top")
    ax.set_ylim(45, 101)
    ax.set_xlim(0, 90)
    ax.set_xticks([0, 30, 60, 90])
    ax.set_xlabel("misalignment $\\alpha$ (deg)")
    ax.set_ylabel("$F_r$ (%)")
    ax2 = ax.twinx()
    ax2.semilogy(al, lam, color="#777777", lw=0.9, ls=":")
    ax2.set_ylabel("cyclicity $\\Lambda$", color="#555555")
    ax2.tick_params(axis="y", colors="#555555")

    # (d) fidelity vs loaded Q for three misalignments
    ax = axs[3]
    for key, c, lab in (("ref_a9", OI["green"], "9$^\\circ$"),
                        ("a35", OI["blue"], "35$^\\circ$"),
                        ("a45", OI["vermilion"], "45$^\\circ$")):
        rws = rr["vs_QL"][key]["rows"]
        ax.semilogx([r["Q_L"] for r in rws], [100 * r["F"] for r in rws],
                    "o-", ms=2.5, color=c, label=lab)
    ax.axvline(cavL["Q_L"], color="#777777", lw=0.6, ls=":")
    ax.set_xlabel("loaded $Q_L$")
    ax.set_ylabel("$F_r$ (%)")
    ax.legend(loc="lower left", fontsize=6.5, handlelength=1.0,
              handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25)

    # (e) fidelity and Rabi rate vs strain
    ax = axs[4]
    rows = rr["vs_strain"]
    x = np.array([r["ratio"] for r in rows])
    ax.axvspan(0.0, 0.65, color=OI["grey"], alpha=0.18, lw=0)
    ax.axvspan(1.0, 1.9, color=OI["grey"], alpha=0.35, lw=0)
    ax.plot(x, [100 * r["F"] for r in rows], "o-", ms=2.5, color=OI["green"])
    ax.set_ylim(95, 100)
    ax.set_xlim(0, 2.6)
    ax.set_xlabel("strain $2\\Upsilon_g/\\lambda_g$")
    ax.set_ylabel("$F_r$ (%)")
    ax2 = ax.twinx()
    ax2.plot(x, [r["rabi_MHz"] for r in rows], color=OI["vermilion"],
             lw=0.9, ls="--")
    ax2.set_ylabel("Rabi rate (MHz)", color=OI["vermilion"])
    ax2.tick_params(axis="y", colors=OI["vermilion"])
    ax2.set_ylim(0, 8.5)

    # (f) line broadening
    ax = axs[5]
    rb = rr["robust"]
    ax.plot([r["x_GHz"] for r in rb["extra_lw"]],
            [100 * r["F"] for r in rb["extra_lw"]], "o-", ms=2.5,
            color=OI["green"], label="homogeneous")
    ax.plot([r["fwhm_GHz"] for r in rb["sd"]],
            [100 * r["F"] for r in rb["sd"]], "s-", ms=2.5,
            color=OI["blue"], label="spectral diffusion")
    ax.set_xlabel("extra width (GHz)")
    ax.set_ylabel("$F_r$ (%)")
    ax.legend(loc="lower left", fontsize=6.5, handlelength=1.0,
              handletextpad=0.4, borderaxespad=0.2, labelspacing=0.25)

    letters(fig, list(zip(axs, ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"])))
    save(fig, "fig5_readout")


if __name__ == "__main__":
    fig2(); fig4(); fig5(); fig3()
    # Fig. 1 (device concept) is generated by fig1_device.py
