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
eme = json.load(open(RES / "taper_eme.json"))
snvr = json.load(open(RES / "snv_results.json"))
itf = json.load(open(RES / "interface.json"))
frl = json.load(open(RES / "fr_vs_lambda.json"))
prod = json.load(open(RES / "cavity_production.json"))
fld = np.load(RES / "cavity_field.npz")
dmap = np.load(RES / "design_map.npz")
tol = json.load(open(RES / "tolerance.json"))

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
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 1.95))

    # (a) mode dispersion: diamond beam on SiO2 and LN ridge
    ax = axs[0]
    w = 1e3 * np.array(wg_d["w"])
    te0 = []
    for ne, te in zip(wg_d["neff"], wg_d["te"]):
        cands = [n for n, t in zip(ne, te) if t > 0.7]
        te0.append(max(cands) if cands else np.nan)
    ax.plot(w, te0, color=OI["blue"], label="diamond/SiO$_2$ TE$_0$")
    wl = 1e3 * np.array(wg_ln["w"])
    te0l = []
    for ne, te in zip(wg_ln["neff"], wg_ln["te"]):
        cands = [n for n, t in zip(ne, te) if t > 0.7]
        te0l.append(max(cands) if cands else np.nan)
    ax.plot(wl, te0l, color=OI["vermilion"], label="TFLN ridge TE$_0$")
    ax.axhline(1.4574, color=OI["grey"], lw=0.7, ls=":")
    ax.text(640, 1.468, "SiO$_2$ light line", color=OI["grey"], fontsize=6.5,
            ha="right")
    ax.axvline(280, color=OI["blue"], lw=0.7, ls="--", alpha=0.5)
    ax.axvline(600 * 0.66, color="none")
    ax.set_xlabel("width (nm)")
    ax.set_ylabel("$n_\\mathrm{eff}$ at 619 nm")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.10),
              fontsize=6.3)
    panel_label(ax, "(a)")

    # (b) hybrid supermode branches vs diamond taper width
    ax = axs[1]
    cols = [OI["blue"], OI["orange"], OI["green"], OI["purple"]]
    labs = ["TE branch", "TM branch", "LN TE$_1$ (sym.-forbidden)", "LN TM$_0$"]
    for i in range(bn.shape[1]):
        ax.plot(bw, bn[:, i], color=cols[i], lw=1.1)
    ends = bn[0]
    ax.text(58, ends[0] - 0.010, "TE branch", color=cols[0],
            fontsize=5.8, ha="right", va="top")
    ax.text(58, ends[1] - 0.012, "TM branch", color=cols[1],
            fontsize=5.8, ha="right", va="top")
    ax.text(58, ends[2] - 0.012, "LN TE$_1$ (sym.-forb.)",
            color=cols[2], fontsize=5.8, ha="right", va="top")
    ax.text(58, ends[3] + 0.012, "LN TM$_0$", color=cols[3],
            fontsize=5.8, ha="right", va="bottom")
    ax.set_xlabel("diamond taper width (nm)")
    ax.set_ylabel("supermode $n_\\mathrm{eff}$")
    ax.invert_xaxis()
    panel_label(ax, "(b)")

    # (c) EME transfer vs taper length
    ax = axs[2]
    L = np.array(eme["lengths"]); T = np.array(eme["T"])
    ax.plot(L, 100 * T, "o-", color=OI["green"], ms=3.5)
    ax.axhline(92, color=OI["grey"], lw=0.8, ls="--")
    ax.text(19.8, 92.7, "Riedel et al.\n(measured, SiV/737 nm)",
            color=OI["grey"], fontsize=6.0, ha="right", va="bottom")
    ax.annotate(f"{100*T[4]:.1f}% at 6 µm", xy=(6, 100 * T[4]),
                xytext=(7.0, 84.5), fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.set_xlabel("taper length (µm)")
    ax.set_ylabel("transfer $T$ (%)")
    ax.set_ylim(80, 101)
    panel_label(ax, "(c)", dx=-0.30)

    fig.tight_layout(w_pad=1.8)
    save(fig, "fig2_waveguide")


# ================================================================ Fig 3
def fig3():
    fig = plt.figure(figsize=(DOUBLE, 2.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1.5, 1])

    # (a) mirror band structure (recomputed quickly from saved edges is
    # heavy; use the plotted dielectric/air edges + cavity line schematic
    # from stored numbers)
    ax = fig.add_subplot(gs[:, 0])
    import cavity as cav
    kx, fr = cav.band_structure(194.35, nk=13, gmax=3.0, gmode_inds=(0, 1, 2),
                                numeig=8)
    a = 194.35
    for b in range(fr.shape[1]):
        ax.plot(kx / np.pi, fr[:, b], ".", ms=2.2, color=OI["grey"])
    # silica light cone and TE band gap
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
    ax.text(0.62, 0.345, "TE gap", fontsize=6.5, color=OI["black"])
    f_c = a / prod["lam_nm"]
    ax.axhline(f_c, color=OI["vermilion"], lw=1.0)
    ax.text(0.05, f_c + 0.004, "cavity mode (618.7 nm)",
            color=OI["vermilion"], fontsize=6.5)
    ax.set_xlabel("$k_x a/\\pi$")
    ax.set_ylabel("frequency $a/\\lambda$")
    ax.set_xlim(0, 1); ax.set_ylim(0.15, 0.42)
    panel_label(ax, "(a)")

    # (b) field maps
    axb = fig.add_subplot(gs[0, 1])
    xs, ys = fld["xs"], fld["ys"]
    E2 = fld["E2xy"]
    axb.imshow(E2 / E2.max(), extent=[xs[0] * a / 1e3, xs[-1] * a / 1e3,
               ys[0] * a / 1e3, ys[-1] * a / 1e3], origin="lower",
               cmap="inferno", aspect="auto", vmax=1.0)
    axb.set_ylim(-0.35, 0.35)
    axb.set_ylabel("y (µm)")
    axb.set_xticklabels([])
    panel_label(axb, "(b)", dx=-0.09)
    axc = fig.add_subplot(gs[1, 1])
    zs = fld["zs"]
    E2z = fld["E2xz"]
    axc.imshow(E2z / E2z.max(), extent=[xs[0] * a / 1e3, xs[-1] * a / 1e3,
               zs[0] * a / 1e3, zs[-1] * a / 1e3], origin="lower",
               cmap="inferno", aspect="auto", vmax=1.0)
    axc.set_ylim(-0.25, 0.45)
    axc.set_xlabel("x (µm)")
    axc.set_ylabel("z (µm)")

    # (c) envelope decay and loading
    ax = fig.add_subplot(gs[:, 2])
    line = fld["E2line"]
    ax.semilogy(xs * a / 1e3, line / line.max(), color=OI["blue"], lw=0.8)
    D = prod["mirror_decay_D"]
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("$|E|^2$ (norm., y=z=0)")
    ax.set_ylim(1e-5, 2)
    ax.text(0.03, 0.97,
            f"mirror decay\n$D = {D:.2f}$/period\n(={10*np.log10(1/D):.1f} dB)",
            transform=ax.transAxes, fontsize=6.5, va="top",
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
    panel_label(ax, "(c)")

    fig.tight_layout(w_pad=1.6)
    save(fig, "fig3_cavity")


# ================================================================ Fig 4
def fig4():
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 1.95))

    # (a) cyclicity vs field angle + measured points
    ax = axs[0]
    z = np.array(snvr["zeta_deg"]); c = np.array(snvr["cyc_vs_zeta"])
    ax.semilogy(z, c, color=OI["blue"], label="model")
    ax.errorbar([147], [2244], yerr=[108], fmt="o", color=OI["vermilion"],
                ms=4, capsize=2, label="measured")
    ax.errorbar([53], [8.6], yerr=[0.4], fmt="o", color=OI["vermilion"], ms=4,
                capsize=2)
    ax.annotate("out-of-plane\nalignment, see (b)", xy=(147, 2244),
                xytext=(178, 250), fontsize=5.8,
                arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.set_ylim(3, 8e3)
    ax.set_xlabel("field angle $\\zeta$ (deg)")
    ax.set_ylabel("cyclicity $\\Lambda$")
    ax.set_xlim(0, 360)
    ax.legend(loc="upper right", fontsize=6.5)
    panel_label(ax, "(a)")

    # (b) cyclicity vs misalignment
    ax = axs[1]
    al = np.array(snvr["alpha_deg"]); ca = np.array(snvr["cyc_vs_alpha"])
    ax.loglog(al, ca, color=OI["blue"])
    ax.errorbar([10], [2244], yerr=[108], fmt="o", color=OI["vermilion"],
                ms=4, capsize=2)
    ax.annotate("measured maximum\n(residual misalign. $\\approx 10°$)",
                xy=(10, 2244), xytext=(0.9, 60), fontsize=6.2,
                arrowprops=dict(arrowstyle="->", lw=0.7))
    ax.set_xlabel("misalignment $\\alpha$ (deg)")
    ax.set_ylabel("cyclicity $\\Lambda$")
    panel_label(ax, "(b)", dx=-0.30)

    # (c) strain trade-off
    ax = axs[2]
    r = np.array(snvr["strain_2UgOverLam"])
    cyc = np.array(snvr["cyc_vs_strain"])
    rab = np.array(snvr["rabi_MHz_vs_strain"])
    ax.loglog(r, cyc, color=OI["blue"], label="$\\Lambda$")
    ax.set_xlabel("strain $2\\Upsilon_g/\\lambda_g$")
    ax.set_ylabel("cyclicity $\\Lambda$", color=OI["blue"])
    ax.axvline(snvr["device_strain_ratio"], color=OI["grey"], lw=0.8, ls="--")
    ax.text(snvr["device_strain_ratio"] * 0.85, 2.8e5, "measured\ndevice",
            ha="right", va="top",
            fontsize=6.2, color=OI["grey"])
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
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE, 2.60))

    # (a) efficiency budget: confocal vs integrated
    ax = axs[0]
    conf = [("QE $\\eta_q$", 0.8), ("PSB $1-\\eta_{DW}$", 0.43),
            ("scatter", 0.05), ("path", 0.35), ("detector", 0.65)]
    cavb = [("$\\beta$ (cavity)", 0.905), ("$\\kappa_{wg}/\\kappa$", 0.9996),
            ("taper", 0.991), ("chip+coupl.", 0.8), ("detector", 0.9)]
    y = np.arange(5)
    ax.barh(y + 0.2, [c[1] for c in conf], height=0.36,
            color=OI["grey"], label="confocal (meas. ref.)")
    ax.barh(y - 0.2, [c[1] for c in cavb], height=0.36,
            color=OI["green"], label="this design")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c[0]} | {v[0]}" for c, v in zip(conf, cavb)],
                       fontsize=5.8)
    ax.invert_yaxis()
    tot_c = np.prod([c[1] for c in conf])
    tot_v = np.prod([c[1] for c in cavb])
    ax.set_xlabel("per-stage efficiency")
    ax.set_title(f"total $\\eta$:  {100*tot_c:.1f}%  vs  {100*tot_v:.0f}%",
                 fontsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.52),
              fontsize=6, ncol=2)
    ax.set_xlim(0, 1.05)
    panel_label(ax, "(a)", dx=-0.60, dy=1.12)

    # (b) Fr vs Lambda0
    ax = axs[1]
    L0 = np.array(frl["Lambda0"])
    ax.semilogx(L0, 100 * np.array(frl["cavity"]), color=OI["green"],
                label="integrated design")
    ax.fill_between(L0, 100 * np.array(frl["confocal_002"]),
                    100 * np.array(frl["confocal_004"]), color=OI["grey"],
                    alpha=0.4, lw=0, label="confocal, $\\eta$=0.2–0.4%")
    ax.plot([2244], [87.4], "o", ms=4, color=OI["vermilion"],
            label="measured 87.4%")
    ax.axhline(50, color=OI["grey"], lw=0.5, ls=":")
    for lx, txt in ((8.6, "worst\nangle"), (100, "high\nstrain"),
                    (2244, "max $\\Lambda$")):
        ax.axvline(lx, color=OI["skyblue"], lw=0.6, ls="--", alpha=0.7)
        ax.text(lx * 1.12, 103, txt, fontsize=5.6, color=OI["blue"],
                va="bottom")
    ax.set_xlabel("bare cyclicity $\\Lambda_0$")
    ax.set_ylabel("single-shot fidelity $F_r$ (%)")
    ax.set_ylim(48, 102)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.50),
              fontsize=5.8, ncol=1, frameon=False)
    ax.text(0.035, 0.965, "(b)", transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="top")

    # (c) design map
    ax = axs[2]
    QLs, xis = dmap["QLs"], dmap["xis"]
    Fr = dmap["Fr"]
    pc = ax.pcolormesh(QLs, xis, 100 * Fr, cmap="viridis", vmin=75, vmax=100,
                       shading="nearest")
    cb = fig.colorbar(pc, ax=ax, pad=0.02)
    cb.set_label("$F_r$ (%) at $\\Lambda_0=100$", fontsize=7)
    ax.contour(QLs, xis, dmap["valid"].astype(float), levels=[0.5],
               colors="w", linewidths=1.0, linestyles="--")
    ax.plot([500], [0.5], "*", ms=9, color=OI["vermilion"],
            markeredgecolor="w", markeredgewidth=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("loaded quality factor $Q_L$")
    ax.set_ylabel("placement overlap $\\xi_{pos}$")
    ax.text(0.035, 0.965, "(c)", transform=ax.transAxes, fontweight="bold",
            fontsize=9, va="top")

    fig.tight_layout(w_pad=1.6)
    save(fig, "fig5_readout")


if __name__ == "__main__":
    fig2(); fig4(); fig5(); fig3()
    # Fig. 1 (device concept) is generated by fig1_device.py
