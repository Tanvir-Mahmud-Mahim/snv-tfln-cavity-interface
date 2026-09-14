"""Stage-4 integration: from photonic design to readout predictions.

Combines the GME cavity (Q_i, V, mirror decay D), the EME taper transfer,
and the validated SnV emission/cyclicity model into the detected-photon
budget and single-shot readout fidelity of the integrated interface, and
maps the design space with its validity boundaries:

  B1 (contrast):        gamma_cav < omega_q / 5   (spin selectivity)
  B2 (bad cavity):      gamma_cav < kappa / 10    (Purcell validity)
  B3 (weak coupling):   g < kappa / 10
  B4 (C-selectivity):   kappa < Delta_CD          (no D enhancement)

Design point: heavily overcoupled one-sided cavity, Q_L ~ 1.5e3,
xi_pos = 0.5, giving F_C ~ 57, beta ~ 0.97, eta_wg ~ 0.999.
"""

from __future__ import annotations

import json
import pathlib

import numpy as np

import readout as ro

RESULTS = pathlib.Path(__file__).resolve().parent / "results"

C0 = 299792458.0
S_SAT = 2.0            # readout drive saturation parameter p/psat

# ---------------------------------------------------- photonic design inputs
LAM_NM = 618.9            # cavity resonance (production GME; tuned to C line)
F_CAV = C0 / (LAM_NM * 1e-9)
Q_I = 1.39e6              # intrinsic (radiation-limited) Q, production GME
V_REL = 0.68              # mode volume, (lam/n)^3, gmax 2.0-3.0 spread
D_MIRROR = 0.399          # per-period intensity decay in the mirror (GME)
T_TAPER = 0.991           # EME taper transfer at L = 6 um
OMEGA_Q = 3.677e9         # qubit frequency (operating point, 2403.13110)
DELTA_CD = 2.097e12       # C/D splitting ~ (Delta_e - Delta_g)
GAMMA_C_FREE = ro.GAMMA_0 * ro.ETA_Q * ro.ETA_DW * ro.ETA_BR / (2 * np.pi)
                          # C-transition free radiative rate /2pi ~ 12.1 MHz

# Fabry-Perot loading model: Q_wg(N) = Q_FP0 / D^N with
# Q_FP0 = 2 pi f tau_rt, tau_rt = 2 n_g L_c / c
N_G, L_C = 10.0, 1.2e-6   # effective group index and cavity length


def q_wg(n_out):
    tau_rt = 2 * N_G * L_C / C0
    return 2 * np.pi * F_CAV * tau_rt / (D_MIRROR ** n_out)


def design(Q_L, xi_pos=0.5, Lambda_0=2244.0, eta_chip=0.8, eta_det=0.9,
           f0=0.99):
    """All derived quantities at a given loaded Q."""
    kappa = F_CAV / Q_L
    F_max = ro.purcell_max(Q_L, V_REL)
    F_C = 1 + (F_max - 1) * ro.XI_POL * xi_pos
    zeta = 1 + ro.ETA_Q * ro.ETA_DW * ro.ETA_BR * (F_C - 1)
    gamma_cav = zeta / ro.TAU_0 / (2 * np.pi)       # linewidth, Hz
    g = 0.5 * np.sqrt(F_C * kappa * GAMMA_C_FREE)   # coupling, Hz
    beta = ro.ETA_Q * ro.ETA_DW * ro.ETA_BR * F_C / zeta
    eta_wg = 1 - Q_L / Q_I
    eta = beta * eta_wg * T_TAPER * eta_chip * eta_det
    # dark-state off-resonant scattering ratio R_dark/R_bright at
    # saturation parameter S_SAT (power broadening included)
    leak = (1.0 + S_SAT) / (1.0 + S_SAT + (2 * OMEGA_Q / gamma_cav) ** 2)
    # D-transition Purcell (cavity Lorentzian at Delta_CD)
    F_D = 1 + (F_max - 1) * ro.XI_POL * xi_pos * ro.lorentzian_suppression(
        DELTA_CD, F_CAV, Q_L)
    return dict(Q_L=Q_L, kappa_GHz=kappa / 1e9, F_max=F_max, F_C=F_C,
                F_D=F_D, zeta=zeta, tau_cav_ns=ro.TAU_0 / zeta * 1e9,
                gamma_cav_GHz=gamma_cav / 1e9, g_GHz=g / 1e9,
                beta=beta, eta_wg=eta_wg, eta=eta, leak=leak,
                ok_contrast=gamma_cav < OMEGA_Q / 5,
                ok_badcav=gamma_cav < kappa / 10,
                ok_weak=g < kappa / 10,
                ok_Csel=kappa < DELTA_CD)


def leak_of(d, s):
    gamma_cav = d["gamma_cav_GHz"] * 1e9
    return (1.0 + s) / (1.0 + s + (2 * OMEGA_Q / gamma_cav) ** 2)


def readout_point(d, Lambda_0, tau=None, noise_rate=1e3, f0=0.99):
    """Readout fidelity at a design point d for bare cyclicity Lambda_0.
    Optimises over the integer threshold, the readout window (1-8
    polarisation times) and the drive saturation parameter s."""
    gamma_cav_rate = d["zeta"] / ro.TAU_0
    best = None
    for s in (0.2, 0.5, 1.0, 2.0, 5.0):
        R = 0.5 * gamma_cav_rate * s / (1 + s)
        Gp = R / (1 + Lambda_0)
        taus = [tau] if tau is not None else list(
            np.array([0.005, 0.02, 0.05, 0.1, 0.3, 1, 2, 3, 5]) / Gp)
        for t in taus:
            Fr, Nr, mb, nd = ro.fidelity_exact(
                d["eta"], Lambda_0, gamma_cav_rate, leak_of(d, s), t,
                noise_rate=noise_rate, s=s, f0=f0)
            if best is None or Fr > best[2]:
                best = (mb, nd, Fr, Gp, t, Nr, s)
    return best


if __name__ == "__main__":
    out = {}
    # loading: mirror periods -> Q_wg
    out["q_wg_vs_N"] = {str(n): float(q_wg(n)) for n in range(0, 11)}

    dp = design(500.0)
    out["design_point"] = {k: (float(v) if not isinstance(v, bool) else v)
                           for k, v in dp.items()}
    print(json.dumps(out["design_point"], indent=1))

    # readout at design point across operating regimes
    regimes = {"max_cyclicity": (2244.0, 2e-6), "high_strain": (100.0, 5e-7),
               "worst_angle": (8.6, 2e-6)}
    out["readout"] = {}
    for name, (lam0, tau) in regimes.items():
        nb, nd, Fr, Gp, topt, Nr, sopt = readout_point(dp, lam0)
        out["readout"][name] = dict(Lambda0=lam0, tau_us=topt * 1e6,
                                    nb=float(nb), nd=float(nd), Fr=float(Fr),
                                    Gp_per_s=float(Gp), Nr=int(Nr),
                                    s=float(sopt))
        print(f"{name:15s} L0={lam0:6.1f} tau={topt*1e6:6.2f}us Nr={Nr:3d} "
              f"s={sopt:.1f} nb={nb:7.1f} nd={nd:.3f} Fr={Fr:.5f}")

    # optimal Q_L per operating cyclicity
    out["QL_opt"] = {}
    for lam0 in (8.6, 100.0, 2244.0):
        scan = []
        for ql in np.logspace(2.3, 4.5, 16):
            d = design(float(ql))
            scan.append((float(ql), readout_point(d, lam0)[2]))
        qlb, frb = max(scan, key=lambda x: x[1])
        out["QL_opt"][str(lam0)] = dict(QL=qlb, Fr=frb)
        print(f"optimal Q_L for Lambda0={lam0}: {qlb:.0f} -> Fr={frb:.5f}")

    # design-space map: Q_L x xi_pos
    import sys
    run_map = "--map" in sys.argv
    QLs = np.logspace(2.3, 5.0, 21)
    xis = np.linspace(0.05, 1.0, 13)
    Fr_map = np.zeros((len(xis), len(QLs)))
    if not run_map:
        QLs_run = []
    else:
        QLs_run = QLs
    eta_map = np.zeros_like(Fr_map)
    valid = np.zeros_like(Fr_map, dtype=bool)
    FC_map = np.zeros_like(Fr_map)
    if run_map:
        for i, xi in enumerate(xis):
            for j, ql in enumerate(QLs):
                d = design(float(ql), float(xi))
                Fr = readout_point(d, 100.0)[2]
                Fr_map[i, j] = Fr
                eta_map[i, j] = d["eta"]
                FC_map[i, j] = d["F_C"]
                valid[i, j] = (d["ok_contrast"] and d["ok_badcav"]
                               and d["ok_weak"])
        np.savez(RESULTS / "design_map.npz", QLs=QLs, xis=xis, Fr=Fr_map,
                 eta=eta_map, FC=FC_map, valid=valid)

    # Fr vs Q_L robustness (fixed xi=0.5, Lambda=100)
    out["Fr_vs_QL"] = {}
    for ql in (500, 1000, 1500, 3000, 5000, 10000):
        d = design(float(ql))
        Fr = readout_point(d, 100.0)[2]
        out["Fr_vs_QL"][str(ql)] = dict(Fr=float(Fr), eta=float(d["eta"]),
                                        F_C=float(d["F_C"]))
    (RESULTS / "interface.json").write_text(json.dumps(out, indent=1))
    print("saved interface.json + design_map.npz")
