"""Production cavity run: save field maps + extract mirror decay constant."""
import json, pathlib
import numpy as np
import legume
import cavity as cav

RES = cav.RESULTS
A = 194.35; GMAX = 2.5; NUMEIG = 130; DEPTH = 0.14; NMIR = 10; WY = 4.0
phc, Lx = cav.cavity_phc(A, N_mirror=NMIR, n_taper=8, taper_depth=DEPTH, W_y=WY)
gme = legume.GuidedModeExp(phc, gmax=GMAX)
gme.run(kpoints=np.array([[0.0],[0.0]]), gmode_inds=[0], numeig=NUMEIG,
        compute_im=False, verbose=False, eps_eff='background')
freqs = np.array(gme.freqs[0])
t_rel = 200.0/A
# identify cavity mode (same criterion as run_cavity)
f_lo = 0.2884  # TE dielectric edge at a=193.4 (from run_cavity output)
in_gap = np.where((freqs > f_lo + 5e-4) & (freqs < f_lo + 0.032))[0]
xs = np.linspace(-Lx/2, Lx/2, 480, endpoint=False)
ys = np.linspace(-WY/2, WY/2, 96, endpoint=False)
best = None
for mi in in_gap:
    fld, _, _ = gme.get_field_xy("e", kind=0, mind=int(mi), z=t_rel/2, xgrid=xs, ygrid=ys)
    E2 = sum(np.abs(fld[c])**2 for c in "xyz")
    frac = E2[:, np.abs(xs) < 3.0].sum() / E2.sum()
    if best is None or frac > best[1]:
        best = (int(mi), float(frac), E2)
mi, frac, E2xy = best
f = float(freqs[mi]); lam = A/f
print(f"cavity mode {mi}: lam={lam:.2f} nm frac={frac:.3f}")
fim = float(np.array(gme.compute_rad(0,[mi])[0]).ravel()[0])
Q = f/(2*fim)
print(f"Q = {Q:.4g}")

# xz map (vertical cross-section along beam axis)
zs = np.linspace(-0.8, t_rel+0.8, 120)
E2xz = np.zeros((len(zs), len(xs)))
for i, z in enumerate(zs):
    fld, _, _ = gme.get_field_xy("e", kind=0, mind=mi, z=z, xgrid=xs, ygrid=np.array([0.0]))
    E2xz[i] = sum(np.abs(fld[c][0])**2 for c in "xyz")

# mirror decay: envelope of E2 at y=0, z=mid, sampled at antinodes between holes
fld, _, _ = gme.get_field_xy("e", kind=0, mind=mi, z=t_rel/2, xgrid=xs, ygrid=np.array([0.0]))
E2line = sum(np.abs(fld[c][0])**2 for c in "xyz")
xs_h, radii = cav.taper_positions(NMIR, 8, DEPTH)
# antinode positions: midpoints between consecutive mirror holes (j>=8)
mids = [(xs_h[j]+xs_h[j+1])/2 for j in range(7, len(xs_h)-1)]
vals = [float(E2line[np.argmin(np.abs(xs - m))]) for m in mids]
decay = [vals[i+1]/vals[i] for i in range(len(vals)-1)]
D = float(np.exp(np.mean(np.log(decay[2:7]))))  # plateau region, away from taper and supercell edge
print("per-period intensity decay factors:", np.round(decay,4), "-> D =", round(D,4))

np.savez(RES/"cavity_field.npz", xs=xs, ys=ys, zs=zs, E2xy=E2xy, E2xz=E2xz,
         E2line=E2line, hole_x=xs_h, hole_r=radii, Lx=Lx, a_nm=A, lam_nm=lam)
out = dict(a_nm=A, gmax=GMAX, lam_nm=lam, Q=Q, mode_index=int(mi),
           center_frac=frac, mirror_decay_D=D, decay_list=[float(d) for d in decay])
(RES/"cavity_production.json").write_text(json.dumps(out, indent=1))
print("saved")
