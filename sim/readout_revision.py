"""Readout predictions of the revised manuscript (Hamiltonian-consistent
four-level model of readout_fcs.py).

Writes results/readout_rev.json.  The loaded-cavity parameters come
from results/cavity_loaded.json (3D FDTD, fdtd_analysis.py) and the taper
transfer from results/tolerance_eme/nominal.json (taper_tolerance.py).
Run:  python readout_revision.py [part ...]   (default: all parts:
reference ql alpha strain field robust extra)
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np
from scipy.optimize import brentq

import snv
import readout_fcs as rf

RES = pathlib.Path(__file__).resolve().parent / "results"
OUT_JSON = RES / "readout_rev.json"

P = snv.PARAMS
TH, PH = P["theta"], P["phi"]
MU = np.array([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)])
EX = np.array([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)])
EY = np.array([-np.sin(PH), np.cos(PH), 0.0])
E1 = np.cross(MU, [0, 0, 1.0]); E1 /= np.linalg.norm(E1)   # = -EY

ALPHA_REF = 9.0            # misalignment reproducing Lambda = 2244 (Sec. 2.3)
FQ_REF = 3.677             # measured qubit frequency (GHz) [Rosenthal 2024]
TAUS = np.geomspace(5.0, 1.0e5, 46)     # 5 ns .. 100 us


def b_dir(alpha_deg, az_deg=None):
    """Unit field direction at polar angle alpha from the SnV axis.
    az_deg=None: the misalignment plane of Fig. 4(b) (towards E1);
    otherwise the azimuth in the spin frame (0 = strain axis x)."""
    a = np.deg2rad(alpha_deg)
    if az_deg is None:
        perp = E1
    else:
        f = np.deg2rad(az_deg)
        perp = np.cos(f) * EX + np.sin(f) * EY
    return np.cos(a) * MU + np.sin(a) * perp


def b_ref():
    """Field magnitude at which the model qubit frequency equals the
    measured 3.677 GHz at the reference misalignment."""
    return brentq(lambda B: snv.qubit_frequency(B * b_dir(ALPHA_REF)) - FQ_REF,
                  0.05, 0.2, xtol=1e-6)


def taper_T():
    """Nominal EME taper transfer (6-um taper, six modes, 30 segments)."""
    d = json.loads((RES / "tolerance_eme" / "nominal.json").read_text())
    return d["T_LN"][d["lengths"].index(6.0)]


def cavity_params():
    """Loaded-cavity parameters: from the 3D FDTD of the one-sided device
    if available (results/cavity_loaded.json), else the analytic model."""
    f = RES / "cavity_loaded.json"
    if f.exists():
        c = json.loads(f.read_text())
        return dict(Q_L=c["Q_L"], V_rel=c["V_eff_equiv"], eta_wg=c["eta_wg"],
                    source="FDTD", table=c["table"])
    return dict(Q_L=500.0, V_rel=0.68, eta_wg=None, source="analytic")


def budget(ls, Q_L=None, xi_pos=0.5, eta_chip=0.5, eta_det=0.77, **kw):
    cp = cavity_params()
    V = cp["V_rel"]
    if Q_L is None:
        Q_L = cp["Q_L"]
        eta_wg = cp["eta_wg"]
    elif "table" in cp:
        # interpolate the FDTD waveguide fraction and equivalent mode
        # volume in log(Q_L) between the simulated N_out points
        t = cp["table"]
        lq = np.log(t["Q_L"])
        eta_wg = float(np.interp(np.log(Q_L), lq, t["eta_wg"]))
        V = float(np.interp(np.log(Q_L), lq, t["V_eff_equiv"]))
    else:
        eta_wg = None
    kw.setdefault("T_taper", taper_T())
    return rf.cavity_budget(Q_L, V_rel=V, xi_pos=xi_pos,
                            delta_cd_ghz=ls["delta_g"], eta_wg=eta_wg,
                            eta_chip=eta_chip, eta_det=eta_det, **kw)


def ls_summary(ls):
    return dict(wq=ls["wq"], we=ls["we"], dA1B2=float(ls["nu"][1, 1]),
                dA1A2=float(ls["nu"][1, 0]), dA1B1=float(ls["nu"][0, 1]),
                LamA=ls["LamA"], LamB=ls["LamB"], delta_g=ls["delta_g"],
                S_rel=(ls["S"] / ls["S"][0, 0]).tolist())


def fid(ls, cb, **kw):
    return rf.readout_fidelity(ls, cb["zeta"] * rf.GAMMA_0, cb["eta"],
                               s_list=(0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
                               taus=TAUS, **kw)


def part_reference(out):
    B = b_ref()
    ls = rf.level_structure(B * b_dir(ALPHA_REF))
    cb = budget(ls)
    r = {"B_T": B, "alpha": ALPHA_REF, "levels": ls_summary(ls),
         "budget": cb, "cavity_source": cavity_params()["source"]}
    best = fid(ls, cb)
    r["best"] = best
    F, so, nr, nb, nd = rf.fidelity_vs_tau(ls, cb["zeta"] * rf.GAMMA_0,
                                           cb["eta"], TAUS)
    r["vs_tau"] = dict(tau_ns=TAUS.tolist(), F=F.tolist(), s=so.tolist(),
                       nb=nb.tolist(), nd=nd.tolist())
    # confocal chain of Rosenthal 2024 at the same field: eta 0.2-0.4%,
    # lifetime-limited and with the measured 60 MHz linewidth of the
    # emitter of Rosenthal 2023 (extra 25 MHz homogeneous), measured
    # background (n_d <= 0.2 counts in 50 us)
    conf = {}
    for eta in (0.002, 0.004):
        for xlw in (0.0, 0.025):
            Fc, sc, _, nbc, ndc = rf.fidelity_vs_tau(
                ls, rf.GAMMA_0, eta, TAUS, gamma_phi_extra=xlw,
                bg_dark=0.2 / 50e3, n_fft=512)
            f50 = rf.fidelity_vs_tau(ls, rf.GAMMA_0, eta, [50e3],
                                     gamma_phi_extra=xlw, bg_dark=0.2 / 50e3,
                                     n_fft=512)
            conf[f"eta{eta}_xlw{xlw}"] = dict(F=Fc.tolist(), F50us=float(f50[0][0]),
                                              nb50=float(f50[3][0]),
                                              nd50=float(f50[4][0]))
    r["confocal"] = conf
    out["reference"] = r


def part_ql(out):
    B = b_ref()
    cases = {"ref_a9": B * b_dir(ALPHA_REF), "a35": B * b_dir(35.0),
             "a45": B * b_dir(45.0),
             "inplane53_180mT": snv.b_lab(53.0, 0.180)}
    QLs = [60, 100, 150, 250, 400, 500, 700, 1000, 1500, 3000]
    r = {}
    for name, Bv in cases.items():
        ls = rf.level_structure(Bv)
        rows = []
        for q in QLs:
            cb = budget(ls, Q_L=float(q))
            b = fid(ls, cb)
            f100 = rf.fidelity_vs_tau(ls, cb["zeta"] * rf.GAMMA_0, cb["eta"],
                                      [100.0])[0][0]
            rows.append(dict(Q_L=q, F=b["F"], tau_ns=b["tau_ns"], s=b["s"],
                             F100ns=float(f100), eta=cb["eta"],
                             zeta=cb["zeta"], F_C=cb["F_C"], F_D=cb["F_D"]))
        r[name] = dict(levels=ls_summary(ls), rows=rows)
        print(name, [(x["Q_L"], round(100 * x["F"], 2)) for x in rows], flush=True)
    out["vs_QL"] = r


def part_alpha(out):
    B = b_ref()
    alphas = [0.5, 2, 5, 9, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 63, 66, 70,
              75, 80, 85, 90]
    rows = []
    for a in alphas:
        ls = rf.level_structure(B * b_dir(a))
        cb = budget(ls)
        b = fid(ls, cb)
        cb2 = budget(ls, Q_L=250.0)
        b2 = fid(ls, cb2)
        rows.append(dict(alpha=a, levels=ls_summary(ls), F=b["F"],
                         tau_ns=b["tau_ns"], F_QL250=b2["F"]))
        print("alpha", a, round(ls["LamA"], 1), round(ls["nu"][1, 1], 3),
              round(100 * b["F"], 2), flush=True)
    out["vs_alpha"] = rows


def part_strain(out):
    B = b_ref()
    ratios = [0.05, 0.1, 0.2, 0.3, 0.43, 0.6, 0.8, 1.0, 1.2, 1.5, 1.9, 2.5]
    rows = []
    for x in ratios:
        scale = x * P["lam_g"] / (2 * P["ups_g"])
        Bv = B * b_dir(ALPHA_REF)
        ls = rf.level_structure(Bv, ups_scale=scale)
        cb = budget(ls)
        b = fid(ls, cb)
        rabi = snv.rabi_rate(Bv, 6e-4 * E1, ups_scale=scale) * 1e3
        rows.append(dict(ratio=x, levels=ls_summary(ls), F=b["F"],
                         tau_ns=b["tau_ns"], rabi_MHz=rabi))
        print("strain", x, round(ls["LamA"], 1), round(ls["nu"][1, 1], 3),
              round(100 * b["F"], 2), round(rabi, 2), flush=True)
    out["vs_strain"] = rows


def part_field(out):
    Bs = [0.04, 0.06, 0.0773, 0.08, 0.10, 0.112, 0.15, 0.2, 0.3, 0.4, 0.6]
    rows = []
    for Bm in Bs:
        ls = rf.level_structure(Bm * b_dir(ALPHA_REF))
        cb = budget(ls)
        b = fid(ls, cb)
        rows.append(dict(B=Bm, levels=ls_summary(ls), F=b["F"],
                         tau_ns=b["tau_ns"]))
        print("B", Bm, round(ls["wq"], 2), round(ls["nu"][1, 1], 3),
              round(100 * b["F"], 2), flush=True)
    out["vs_field"] = rows


def part_robust(out):
    B = b_ref()
    ls = rf.level_structure(B * b_dir(ALPHA_REF))
    cb = budget(ls)
    g = cb["zeta"] * rf.GAMMA_0
    r = {}
    r["extra_lw"] = []
    for x in (0.0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5):
        r["extra_lw"].append(dict(x_GHz=x, F=fid(ls, cb, gamma_phi_extra=x)["F"]))
        print("xlw", x, r["extra_lw"][-1]["F"], flush=True)
    r["sd"] = []
    for x in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5):
        r["sd"].append(dict(fwhm_GHz=x, F=fid(ls, cb, sd_fwhm=x)["F"]))
        print("sd", x, r["sd"][-1]["F"], flush=True)
    r["sbr"] = []
    for x in (3, 10, 30, 90, 300, 1000, np.inf):
        r["sbr"].append(dict(sbr=float(x) if np.isfinite(x) else None,
                             F=fid(ls, cb, sbr=x)["F"]))
        print("sbr", x, r["sbr"][-1]["F"], flush=True)
    r["eta_chain"] = []
    for prod in (0.1, 0.2, 0.3, 0.385, 0.5, 0.7, 0.82):
        cbx = budget(ls, eta_chip=prod, eta_det=1.0)
        r["eta_chain"].append(dict(eta_chip_det=prod, eta=cbx["eta"],
                                   F=fid(ls, cbx)["F"]))
        print("eta", prod, r["eta_chain"][-1]["F"], flush=True)
    # combined realistic case
    comb = fid(ls, cb, gamma_phi_extra=0.025, sd_fwhm=0.1, sbr=90.0)
    r["combined"] = comb
    print("combined", comb, flush=True)
    out["robust"] = r


def part_extra(out):
    B = b_ref()
    ls = rf.level_structure(B * b_dir(ALPHA_REF))
    cb = budget(ls)
    g = cb["zeta"] * rf.GAMMA_0
    r = {}
    # placement overlap: best F over Q_L for each xi_pos
    r["xi_pos"] = []
    for xi in (0.1, 0.25, 0.5, 1.0):
        best = None
        for q in (60, 100, 150, 250, 400, 500, 700, 1000, 1500, 3000):
            c = budget(ls, Q_L=float(q), xi_pos=xi)
            f = fid(ls, c)["F"]
            if best is None or f > best[1]:
                best = (q, f)
        r["xi_pos"].append(dict(xi_pos=xi, QL_opt=best[0], F=best[1]))
        print("xi", xi, best, flush=True)
    # D-line-only detection (spectrally separated from the laser)
    cp = cavity_params()
    num_D = (1 - rf.ETA_BR) * (cb["F_D"] - 1)
    num_C = rf.ETA_BR * (cb["F_C"] - 1)
    fD = num_D / (num_C + num_D)
    etaD = cb["eta"] * fD
    r["dline"] = dict(frac=fD, eta=etaD,
                      F=rf.readout_fidelity(ls, g, etaD,
                                            s_list=(0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
                                            taus=TAUS)["F"])
    # direct diamond-to-fiber coupling (91 %, Burek 2017) replaces taper
    # transfer and LN-to-fiber coupling
    etaF = cb["beta"] * cb["eta_wg"] * 0.91 * 0.77     # no taper, no LN coupler
    r["fiber"] = dict(eta=etaF, F=fid(ls, dict(cb, eta=etaF))["F"])
    # A1-B2 detuning equal to the reported 500 MHz
    B500 = brentq(lambda Bm: rf.level_structure(Bm * b_dir(ALPHA_REF))["nu"][1, 1] - 0.5,
                  0.03, 0.12)
    ls5 = rf.level_structure(B500 * b_dir(ALPHA_REF))
    r["B500"] = dict(B=B500, F=fid(ls5, budget(ls5))["F"])
    # conditional rates at the optimal drive of the 100-ns window
    F100, s100, _, nb100, nd100 = rf.fidelity_vs_tau(ls, g, cb["eta"], [100.0])
    cr = rf.conditional_rates(ls, g, float(s100[0]))
    r["cond_rates"] = dict(s=float(s100[0]), F100=float(F100[0]),
                           nb100=float(nb100[0]), nd100=float(nd100[0]),
                           **{k: float(v) for k, v in cr.items()})
    # azimuthal invariance check (V16)
    Fs = []
    for az in (0.0, 45.0, 90.0, 135.0):
        lsa = rf.level_structure(B * b_dir(30.0, az))
        Fs.append(fid(lsa, budget(lsa))["F"])
    r["azimuth_F"] = Fs
    print(json.dumps(r, indent=1, default=float), flush=True)
    out["extra"] = r


PARTS = dict(reference=part_reference, ql=part_ql, alpha=part_alpha,
             strain=part_strain, field=part_field, robust=part_robust,
             extra=part_extra)

if __name__ == "__main__":
    out = json.loads(OUT_JSON.read_text()) if OUT_JSON.exists() else {}
    parts = sys.argv[1:] or list(PARTS)
    for p in parts:
        t = time.time()
        PARTS[p](out)
        out.setdefault("_runtime_s", {})[p] = time.time() - t
        OUT_JSON.write_text(json.dumps(out, indent=1, default=float))
        print(f"[{p}] done in {time.time()-t:.0f}s", flush=True)
