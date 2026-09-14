"""Diamond nanobeam photonic-crystal cavity by guided-mode expansion.

Solver: legume (guided-mode expansion, GME), Minkov et al.,
"Inverse design of photonic crystals through automatic differentiation,"
ACS Photonics 7, 1729 (2020) -- open source, https://github.com/fancompute/legume

Geometry: a one-dimensional photonic-crystal nanobeam etched into a
single-crystal diamond membrane of thickness t = 200 nm and width
w = 280 nm (the single-TE-mode waveguide of ``waveguide.py``), resting on
SiO2 (lower cladding) with air above.  Circular holes of radius r = r_rel*a.
The cavity is formed by a deterministic quadratic taper of the lattice
constant over ``n_taper`` cells on each side of the centre, from the mirror
period ``a`` down to ``a_min = (1-taper_depth)*a``, following the
deterministic high-Q design approach of Quan & Loncar, Opt. Express 19,
18529 (2011), doi:10.1364/OE.19.018529.

The nanobeam is modelled in a supercell: periodic in x with the full
device length, and a wide artificial period W_y in y (air on both sides of
the beam).  All lengths inside legume are in units of the mirror period a;
frequencies come out in units c/a.

Outputs (written to results/):
  cavity_bands.json   mirror-cell band structure and band gap
  cavity_mode.json    cavity frequency, Q, V_eff, Purcell factor
  cavity_conv.json    convergence versus gmax / W_y / N_mirror
  cavity_field.npz    |E|^2 field maps for figures
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import legume

import materials

RESULTS = pathlib.Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

LAM = materials.LAMBDA_SNV_UM                # 0.619 um
N_DIA = float(materials.diamond.n(LAM))
N_SIO2 = float(materials.silica.n(LAM))
EPS_DIA = N_DIA ** 2
EPS_SIO2 = N_SIO2 ** 2

# --- design parameters (lengths in units of the mirror period a) ---------
R_REL = 0.30      # hole radius r = 0.30 a

# legume superposes the Fourier transforms of overlapping shapes as
# sum_i (eps_i - eps_b) * mask_i, so overlapping shapes must be avoided.
# The layer background is therefore set to diamond, and air is added as
# non-overlapping shapes: two rectangles flanking the beam plus the holes.


def _add_air_sides(layer, half_x: float, w_rel: float, W_y: float):
    """Air rectangles flanking the diamond beam (background = diamond)."""
    for lo, hi in ((w_rel / 2, W_y / 2), (-W_y / 2, -w_rel / 2)):
        layer.add_shape(
            legume.Poly(eps=1.0,
                        x_edges=[-half_x, half_x, half_x, -half_x],
                        y_edges=[lo, lo, hi, hi]))


def mirror_cell(a_nm: float, r_rel: float = R_REL, w_nm: float = 280.0,
                t_nm: float = 200.0, W_y: float = 4.0, suspended=False):
    """Single mirror unit cell of the nanobeam (period 1 in x)."""
    lattice = legume.Lattice([1.0, 0], [0, W_y])
    eps_l = 1.0 if suspended else EPS_SIO2
    phc = legume.PhotCryst(lattice, eps_l=eps_l, eps_u=1.0)
    phc.add_layer(d=t_nm / a_nm, eps_b=EPS_DIA)
    _add_air_sides(phc.layers[-1], 0.5, w_nm / a_nm, W_y)
    phc.layers[-1].add_shape(legume.Circle(eps=1.0, x_cent=0, y_cent=0,
                                           r=r_rel))
    return phc


def band_structure(a_nm: float, nk: int = 13, gmax: float = 4.0,
                   gmode_inds=(0, 1, 2), numeig: int = 8,
                   eps_eff="background", **kw):
    phc = mirror_cell(a_nm, **kw)
    gme = legume.GuidedModeExp(phc, gmax=gmax)
    kx = np.linspace(0.0, np.pi, nk)
    kpoints = np.vstack((kx, np.zeros_like(kx)))
    gme.run(kpoints=kpoints, gmode_inds=list(gmode_inds), numeig=numeig,
            compute_im=False, verbose=False, eps_eff=eps_eff)
    return kx, np.array(gme.freqs)   # freqs in c/a


def taper_positions(N_mirror: int, n_taper: int, taper_depth: float):
    """Hole x-positions (units of mirror period a) for a symmetric cavity.

    Quadratic taper of the local period: cell j (j=0 nearest centre) has
    period a_j = a*(1 - taper_depth*(1 - (j/n_taper)^2)) for j < n_taper,
    then N_mirror mirror cells with period a.  Hole radii scale with the
    local period (r_j = r_rel * a_j) to keep r/a fixed.
    """
    periods = []
    for j in range(n_taper):
        periods.append(1.0 - taper_depth * (1.0 - (j / n_taper) ** 2))
    periods += [1.0] * N_mirror
    xs, radii = [], []
    x = periods[0] / 2.0          # first hole half a (tapered) period out
    for j, p in enumerate(periods):
        if j > 0:
            x += (periods[j - 1] + p) / 2.0
        xs.append(x)
        radii.append(R_REL * p)
    return np.array(xs), np.array(radii)


def cavity_phc(a_nm: float, N_mirror: int = 10, n_taper: int = 8,
               taper_depth: float = 0.14, w_nm: float = 280.0,
               t_nm: float = 200.0, W_y: float = 4.0, pad: float = 0.0,
               suspended=False):
    xs, radii = taper_positions(N_mirror, n_taper, taper_depth)
    Lx = 2 * (xs[-1] + 0.5) + 2 * pad
    lattice = legume.Lattice([Lx, 0], [0, W_y])
    eps_l = 1.0 if suspended else EPS_SIO2
    phc = legume.PhotCryst(lattice, eps_l=eps_l, eps_u=1.0)
    phc.add_layer(d=t_nm / a_nm, eps_b=EPS_DIA)
    _add_air_sides(phc.layers[-1], Lx / 2, w_nm / a_nm, W_y)
    for x, r in zip(xs, radii):
        for sgn in (+1, -1):
            phc.layers[-1].add_shape(
                legume.Circle(eps=1.0, x_cent=sgn * x, y_cent=0.0, r=r))
    return phc, Lx


def run_cavity(a_nm: float, gmax: float = 2.0, gmode_inds=(0,),
               numeig: int = 4, f_target: float = None, **kw):
    """Solve the cavity supercell at k=0; return (gme, phc, Lx)."""
    phc, Lx = cavity_phc(a_nm, **kw)
    gme = legume.GuidedModeExp(phc, gmax=gmax)
    opts = dict(gmode_inds=list(gmode_inds), numeig=numeig,
                compute_im=True, verbose=False)
    if f_target is not None:
        opts.update(eig_solver="eigsh", eig_sigma=f_target,
                    use_sparse=False)
    gme.run(kpoints=np.array([[0.0], [0.0]]), **opts)
    return gme, phc, Lx


def mode_metrics(gme, phc, a_nm: float, mode_ind: int, w_nm=280.0,
                 t_nm=200.0, Lx=None, W_y=4.0, nz: int = 21, fim=None):
    """Q, V_eff (in um^3 and (lam/n)^3), Purcell factor."""
    f = float(gme.freqs[0, mode_ind])            # c/a units
    if fim is None:
        fim = float(gme.freqs_im[0, mode_ind])
    Q = f / (2 * fim) if fim > 0 else np.inf
    lam_nm = a_nm / f
    # 3D |E|^2 * eps integral on a grid
    t_rel = t_nm / a_nm
    zs = np.linspace(1e-3, t_rel - 1e-3, nz)     # inside the slab layer
    # legume z-coordinate: layer occupies z in [0, d] (claddings outside)
    Nx, Ny = 220, 60
    xs = np.linspace(-Lx / 2, Lx / 2, Nx, endpoint=False)
    ys = np.linspace(-W_y / 2, W_y / 2, Ny, endpoint=False)
    ueps_max, integ = 0.0, 0.0
    dz = zs[1] - zs[0]
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]
    for z in zs:
        fields, _, _ = gme.get_field_xy("e", kind=0, mind=mode_ind, z=z,
                                        xgrid=xs, ygrid=ys)
        E2 = (np.abs(fields["x"]) ** 2 + np.abs(fields["y"]) ** 2 +
              np.abs(fields["z"]) ** 2)
        eps_r = np.array(phc.get_eps((xs[None, :] * np.ones((Ny, 1)),
                                      ys[:, None] * np.ones((1, Nx)),
                                      z * np.ones((Ny, Nx)))))
        u = eps_r * E2
        integ += float(np.sum(u)) * dx * dy * dz
        ueps_max = max(ueps_max, float(np.max(u)))
    V_a3 = integ / ueps_max                       # in units a^3
    V_um3 = V_a3 * (a_nm * 1e-3) ** 3
    lam_n3 = (lam_nm * 1e-3 / N_DIA) ** 3         # (lam/n)^3 in um^3
    V_rel = V_um3 / lam_n3
    Fp = 3.0 / (4 * np.pi ** 2) * Q / V_rel
    return dict(f_ca=f, lam_nm=lam_nm, Q=Q, V_um3=V_um3,
                V_lam_n3=V_rel, Purcell=Fp)


if __name__ == "__main__":
    # quick mirror-cell band structure scan to centre the gap at 619 nm
    for a in (170, 180, 190, 200):
        kx, fr = band_structure(a, nk=7, gmax=3.0)
        # TE-like fundamental bands at band edge kx=pi
        edge = fr[-1, :4]
        print(f"a={a} nm  band-edge f(c/a)={np.round(edge, 4)}  "
              f"-> lam_edge={np.round(a/edge,1)} nm", flush=True)
