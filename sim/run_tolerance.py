"""Fabrication tolerance + suspended comparison, sequential GME runs."""
import json, numpy as np, legume
import cavity as cav
RES = cav.RESULTS
def solve_one(tag, a_nm=194.35, r_rel=0.30, w_nm=280.0, t_nm=200.0, susp=False):
    cav.R_REL, old = r_rel, cav.R_REL
    phc, Lx = cav.cavity_phc(a_nm, N_mirror=10, n_taper=8, taper_depth=0.14,
                             w_nm=w_nm, t_nm=t_nm, W_y=4.0, suspended=susp)
    cav.R_REL = old
    gme = legume.GuidedModeExp(phc, gmax=2.0)
    gme.run(kpoints=np.array([[0.0],[0.0]]), gmode_inds=[0], numeig=130,
            compute_im=False, verbose=False, eps_eff='background')
    freqs = np.array(gme.freqs[0]); t_rel = t_nm/a_nm
    xs = np.linspace(-Lx/2, Lx/2, 240, endpoint=False)
    ys = np.linspace(-2.0, 2.0, 40, endpoint=False)
    best = None
    for mi in range(40, 130):
        f = freqs[mi]
        if not (0.28 < f < 0.325): continue
        fld,_,_ = gme.get_field_xy("e", kind=0, mind=mi, z=t_rel/2, xgrid=xs, ygrid=ys)
        E2 = sum(np.abs(fld[c])**2 for c in "xyz")
        frac = E2[:, np.abs(xs)<3.0].sum()/E2.sum()
        if best is None or frac > best[1]: best = (mi, frac)
    mi, frac = best
    fim = float(np.array(gme.compute_rad(0,[mi])[0]).ravel()[0])
    out = dict(tag=tag, a_nm=a_nm, r_rel=r_rel, w_nm=w_nm, t_nm=t_nm,
               suspended=susp, lam_nm=float(a_nm/freqs[mi]),
               Q=float(freqs[mi]/(2*fim)), frac=float(frac), mode=int(mi))
    print(json.dumps(out), flush=True)
    return out
res = []
res.append(solve_one("nominal_g2"))
res.append(solve_one("r_plus", r_rel=0.31))
res.append(solve_one("r_minus", r_rel=0.29))
res.append(solve_one("w_plus", w_nm=290.0))
res.append(solve_one("w_minus", w_nm=270.0))
res.append(solve_one("t_minus", t_nm=190.0))
res.append(solve_one("suspended", susp=True))
json.dump(res, open(RES/"tolerance.json","w"), indent=1)
