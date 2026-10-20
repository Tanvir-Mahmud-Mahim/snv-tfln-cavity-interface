"""Analyse the loaded-cavity FDTD runs (cavity_fdtd.py) and write
results/cavity_loaded.json (used by readout_revision.py) and
results/cavity_fdtd_summary.json.

For each N_out:
  * resonance f0 and loaded Q_L from Harminv on the ring-down;
  * Purcell factor of the cavity mode, F_cav = 3 Q_L (lambda/n)^3 /
    (4 pi^2 V), with V the mode volume at the dipole sampled during the
    ring-down (valid also when the ring-down is not complete);
  * fractions of the emitted power at f0 coupled into the TE0 mode of
    the output waveguide, through the closed-side mirror, upward (air),
    downward (substrate) and laterally (ratios of the same ringing field);
  * the non-resonant LDOS background F_bg (relative to bulk diamond)
    away from the resonance, and, where the ring-down is complete
    (N_out = 0), a Lorentzian fit of the LDOS as an independent check of
    F_cav and Q_L;
  * eta_wg, the fraction of the CAVITY-MODE decay that reaches the
    waveguide TE0 mode: frac_TE0 / [F_cav / (F_cav + F_bg)].
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
from scipy.optimize import curve_fit

RES = pathlib.Path(__file__).resolve().parent / "results"
N_DIA = 2.4137


def lor(f, Fbg, Fc, f0, Q):
    return Fbg + Fc / (1 + (2 * Q * (f - f0) / f0) ** 2)


def analyse(n_out, res=40):
    d = json.loads((RES / f"cavity_fdtd_N{n_out}_r{res}.json").read_text())
    f = np.array(d["freqs"])                       # flux grid (101 points)
    fl = {k: np.array(v) for k, v in d["flux"].items()}
    tot = sum(fl.values())
    hm = sorted(d["harminv"], key=lambda m: -m[2])
    f0, Q = hm[0][0], hm[0][1]
    lam = 1.0 / f0
    # Purcell factor of the cavity mode for an ideally oriented dipole at
    # the antinode, from the mode volume sampled during the ring-down
    V_rel = d["V_um3"] / (lam / N_DIA) ** 3
    F_cav = 3.0 / (4.0 * np.pi ** 2) * Q / V_rel
    out = dict(n_out=n_out, resolution=res, f0=f0, lambda_nm=1e3 * lam,
               Q_L=Q, V_rel=V_rel, F_cav=F_cav, t_after=d.get("t_after"))
    # power fractions at f0 (ratios of the same ringing field, robust to
    # an incomplete ring-down)
    at = lambda arr: float(np.interp(f0, f, arr))
    Ptot = at(tot)
    out.update(frac_TE0=at(np.array(d["P_te0"])) / Ptot,
               frac_right_total=at(fl["right"]) / Ptot,
               frac_left=at(fl["left"]) / Ptot, frac_up=at(fl["up"]) / Ptot,
               frac_down=at(fl["down"]) / Ptot,
               frac_side=(at(fl["yp"]) + at(fl["ym"])) / Ptot)
    # non-resonant LDOS background (far from the resonance) and, where
    # the ring-down is complete, the Lorentzian fit of the LDOS spectrum
    bfile = RES / f"cavity_fdtd_bulk_r{res}.json"
    if not bfile.exists():
        # no bulk run at this resolution: normalize with the resolution-40
        # bulk LDOS (the bulk LDOS of a point dipole changes little with
        # resolution; used only for the small non-resonant background)
        bfile = RES / "cavity_fdtd_bulk_r40.json"
        out["bulk_ref_res"] = 40
    if bfile.exists():
        b = json.loads(bfile.read_text())
        ld = np.array(d["ldos"])
        fL = np.linspace(f.min(), f.max(), len(ld))
        F = ld / np.interp(fL, np.array(b["freqs"]), np.array(b["ldos"]))
        far = np.abs(fL - f0) > 10 * f0 / Q
        F_bg = float(np.median(F[far]))
        out["F_bg"] = F_bg
        out["F_ldos_peak"] = float(F.max())
        try:
            popt, _ = curve_fit(lor, fL, F, p0=[F_bg, F.max(), f0, Q],
                                maxfev=20000)
            out.update(F_cav_ldos=float(popt[1]), Q_ldos=float(popt[3]),
                       F_bg_fit=float(popt[0]))
        except Exception:
            pass
    else:
        raise FileNotFoundError("bulk LDOS reference cavity_fdtd_bulk_r40.json missing")
    Fb = out["F_bg"]
    out["cav_frac"] = F_cav / (F_cav + Fb)
    out["eta_wg"] = out["frac_TE0"] / out["cav_frac"]
    return out


if __name__ == "__main__":
    summ = {}
    for n in range(0, 6):
        for res in (40, 50):
            p = RES / f"cavity_fdtd_N{n}_r{res}.json"
            if p.exists():
                a = analyse(n, res)
                summ[f"N{n}_r{res}"] = a
                print(f"N_out={n} res={res}:", json.dumps(a))
    (RES / "cavity_fdtd_summary.json").write_text(json.dumps(summ, indent=1))
    # design point: the most strongly overcoupled termination simulated
    # (smallest N_out); the readout optimum in Q_L is flat at low loaded Q
    r40 = {int(k[1]): v for k, v in summ.items() if k.endswith("_r40")}
    if r40:
        Ns = sorted(r40)
        nd = min(Ns)
        d = r40[nd]
        # equivalent volume such that 3 Q/(4 pi^2 V) - 1 = F_cav, the form
        # used by readout_fcs.cavity_budget (F_max - 1 = F_cav)
        veq = lambda a: 3.0 * a["Q_L"] / (4 * np.pi ** 2 * (a["F_cav"] + 1.0))
        tap = json.loads((RES / "tolerance_eme" / "nominal.json").read_text())
        out = dict(n_out=nd, Q_L=d["Q_L"], F_cav=d["F_cav"],
                   eta_wg=d["eta_wg"], V_eff_equiv=veq(d),
                   lambda_nm=d["lambda_nm"],
                   T_taper=tap["T_LN"][tap["lengths"].index(6.0)],
                   table=dict(n_out=Ns, Q_L=[r40[n]["Q_L"] for n in Ns],
                              eta_wg=[r40[n]["eta_wg"] for n in Ns],
                              V_eff_equiv=[veq(r40[n]) for n in Ns]))
        (RES / "cavity_loaded.json").write_text(json.dumps(out, indent=1))
        print("design:", json.dumps(out))
