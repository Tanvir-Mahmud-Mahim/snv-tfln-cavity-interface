"""Solve the nanobeam cavity, identify the fundamental mode, compute Q/V/Purcell."""
import json, pathlib, sys
import numpy as np
import legume
import cavity as cav

RES = cav.RESULTS
A_NM = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
GMAX = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
NUMEIG = int(sys.argv[3]) if len(sys.argv) > 3 else 130
DEPTH = float(sys.argv[4]) if len(sys.argv) > 4 else 0.14
NMIR = int(sys.argv[5]) if len(sys.argv) > 5 else 10
WY = float(sys.argv[6]) if len(sys.argv) > 6 else 4.0

phc, Lx = cav.cavity_phc(A_NM, N_mirror=NMIR, n_taper=8, taper_depth=DEPTH, W_y=WY)
print(f"a={A_NM} nm Lx={Lx:.2f} a, gmax={GMAX}, numeig={NUMEIG}", flush=True)
gme = legume.GuidedModeExp(phc, gmax=GMAX)
gme.run(kpoints=np.array([[0.0],[0.0]]), gmode_inds=[0], numeig=NUMEIG,
        compute_im=False, verbose=False, eps_eff='background')
freqs = np.array(gme.freqs[0])
print("top 12 freqs:", np.round(freqs[-12:], 4), flush=True)

# gap window from mirror cell (recompute quickly)
kx, fr = cav.band_structure(A_NM, nk=3, gmax=4.0, gmode_inds=(0,1,2), numeig=8)
f_lo = float(fr[-1, 0])        # TE dielectric band edge (Ey-dominant, verified)
print(f"TE dielectric edge: {f_lo:.4f} (c/a) -> {A_NM/f_lo:.1f} nm", flush=True)
# lattice-taper cavity modes are pulled up from the dielectric edge:
in_gap = np.where((freqs > f_lo + 5e-4) & (freqs < f_lo + 0.025))[0]
print("modes in gap:", in_gap, np.round(freqs[in_gap], 4), flush=True)
if len(in_gap) == 0:
    print("NO GAP MODES - increase numeig"); sys.exit(1)

# localization check: |E|^2 at slab mid-plane, energy inside |x|<3 vs total
t_rel = 200.0 / A_NM
xs = np.linspace(-Lx/2, Lx/2, 240, endpoint=False)
ys = np.linspace(-WY/2, WY/2, 40, endpoint=False)
best = None
for mi in in_gap:
    fields, _, _ = gme.get_field_xy("e", kind=0, mind=int(mi), z=t_rel/2,
                                    xgrid=xs, ygrid=ys)
    E2 = sum(np.abs(fields[c])**2 for c in "xyz")
    frac = E2[:, np.abs(xs) < 3.0].sum() / E2.sum()
    print(f"  mode {mi}: f={freqs[mi]:.4f} lam={A_NM/freqs[mi]:.1f} nm center-frac={frac:.3f}", flush=True)
    if best is None or frac > best[1]:
        best = (int(mi), float(frac))
mi, frac = best
print(f"cavity mode = {mi} (center fraction {frac:.3f})", flush=True)

fim = gme.compute_rad(0, [mi])[0]
fim = float(np.array(fim).ravel()[0])
m = cav.mode_metrics(gme, phc, A_NM, mi, Lx=Lx, W_y=WY, fim=fim)
m.update(a_nm=A_NM, gmax=GMAX, numeig=NUMEIG, taper_depth=DEPTH,
         N_mirror=NMIR, W_y=WY, mode_index=mi, center_frac=frac,
         gap_lo_ca=f_lo)
print(json.dumps(m, indent=1), flush=True)
tag = f"a{A_NM:.0f}_g{GMAX}_d{DEPTH}_N{NMIR}_W{WY}"
(RES / f"cavity_mode_{tag}.json").write_text(json.dumps(m, indent=1))
