"""3D FDTD of the LOADED one-sided nanobeam cavity (revision 1,
Reviewer 2 point 6).

The GME calculation of cavity.py uses a closed, symmetric supercell and
therefore gives only the intrinsic (radiation-limited) Q.  Here the
actual one-sided device is simulated with Meep (finite-difference time
domain, open boundaries): ten mirror periods plus the eight-cell
quadratic taper on the closed side, the eight-cell taper plus N_out
mirror periods on the output side, after which the unpatterned
280 x 200 nm diamond beam continues as the output waveguide into a
perfectly matched layer.  The beam rests on a SiO2 half-space with air
above.  A y-polarised point dipole at the field antinode (beam centre)
excites the mode.

Outputs (results/cavity_fdtd_N{N}_r{res}.json): frequencies, the flux
spectra through the six faces of a box enclosing the device, the
eigenmode coefficient of the TE0 mode of the output waveguide (MPB), the
dipole LDOS spectrum, the Harminv modes (f, Q, amplitude), and the mode
volume V = int eps|E|^2 / (eps(r0)|E_y(r0)|^2) sampled from the ringing
field.  fdtd_analysis.py turns these into Q_L, the Purcell factor and
the power fractions.  "bulk" runs the same dipole in bulk diamond (LDOS
reference).

Usage:  python cavity_fdtd.py N_out|bulk [resolution] [ring-down time]
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np
import meep as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from cavity import taper_positions  # noqa: E402
import materials  # noqa: E402

RES_DIR = pathlib.Path(__file__).resolve().parent / "results"
RES_DIR.mkdir(exist_ok=True)
LAM = materials.LAMBDA_SNV_UM
N_DIA = float(materials.diamond.n(LAM))
N_SIO2 = float(materials.silica.n(LAM))

A = 0.1943          # mirror period (um), production design
R_REL = 0.30
W = 0.280
T = 0.200
N_TAPER, DEPTH, N_MIRROR = 8, 0.14, 10


def geometry(n_out, x_right_end):
    dia = mp.Medium(index=N_DIA)
    sio2 = mp.Medium(index=N_SIO2)
    geo = [mp.Block(center=mp.Vector3(0, 0, -5.0), size=mp.Vector3(mp.inf, mp.inf, 10.0),
                    material=sio2),
           mp.Block(center=mp.Vector3(0, 0, T / 2), size=mp.Vector3(mp.inf, W, T),
                    material=dia)]
    xsL, rL = taper_positions(N_MIRROR, N_TAPER, DEPTH)
    xsR, rR = taper_positions(n_out, N_TAPER, DEPTH)
    for x, r in zip(xsL, rL):
        geo.append(mp.Cylinder(radius=r * A, height=T, center=mp.Vector3(-x * A, 0, T / 2),
                               material=mp.air))
    for x, r in zip(xsR, rR):
        geo.append(mp.Cylinder(radius=r * A, height=T, center=mp.Vector3(x * A, 0, T / 2),
                               material=mp.air))
    return geo, xsL[-1] * A, xsR[-1] * A


def run(n_out, resolution=40, dpml=0.5, pad_y=0.45, pad_up=0.45, pad_dn=0.45,
        wg_len=1.0, f_guess=None, df_rel=0.05, t_after=450.0):
    t0 = time.time()
    xL, xR = taper_positions(N_MIRROR, N_TAPER, DEPTH)[0][-1] * A, \
        taper_positions(n_out, N_TAPER, DEPTH)[0][-1] * A
    x_min = -(xL + 0.25) - dpml
    x_max = xR + wg_len + dpml
    sx = x_max - x_min
    cx = 0.5 * (x_max + x_min)
    sy = W + 2 * pad_y + 2 * dpml
    z_min = -pad_dn - dpml
    z_max = T + pad_up + dpml
    sz = z_max - z_min
    cz = 0.5 * (z_max + z_min)
    geo, _, _ = geometry(n_out, x_max)
    # shift the coordinate frame so the cell is centred
    shift = mp.Vector3(-cx, 0, -cz)
    for g in geo:
        g.center = g.center + shift
    src_pt = mp.Vector3(0, 0, T / 2) + shift

    fcen = f_guess if f_guess is not None else 1.0 / LAM
    df = df_rel * fcen
    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey,
                         center=src_pt)]
    sym = [mp.Mirror(mp.Y, phase=-1)]
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy, sz), geometry=geo,
                        sources=sources, resolution=resolution,
                        boundary_layers=[mp.PML(dpml)], symmetries=sym,
                        eps_averaging=True)

    # flux monitors (planes just inside the PML)
    xin_L = x_min + dpml + 0.05
    xin_R = x_max - dpml - 0.05
    yin = sy / 2 - dpml - 0.05
    zt = z_max - dpml - 0.05
    zb = z_min + dpml + 0.05
    nf = 101
    def fr(center, size, weight=1.0):
        return mp.FluxRegion(center=mp.Vector3(*center) + shift, size=mp.Vector3(*size),
                             weight=weight)
    lx = xin_R - xin_L
    ly = 2 * yin
    lz = zt - zb
    xc = 0.5 * (xin_L + xin_R)
    zc = 0.5 * (zt + zb)
    f_right = sim.add_flux(fcen, df, nf, fr((xin_R, 0, zc), (0, ly, lz)))
    f_left = sim.add_flux(fcen, df, nf, fr((xin_L, 0, zc), (0, ly, lz), -1.0))
    f_up = sim.add_flux(fcen, df, nf, fr((xc, 0, zt), (lx, ly, 0)))
    f_dn = sim.add_flux(fcen, df, nf, fr((xc, 0, zb), (lx, ly, 0), -1.0))
    f_yp = sim.add_flux(fcen, df, nf, fr((xc, yin, zc), (lx, 0, lz)))
    f_ym = sim.add_flux(fcen, df, nf, fr((xc, -yin, zc), (lx, 0, lz), -1.0))
    # mode monitor on the output waveguide (smaller window around the beam)
    mode_reg = mp.FluxRegion(center=mp.Vector3(xin_R, 0, T / 2) + shift,
                             size=mp.Vector3(0, ly, lz))
    f_mode = sim.add_mode_monitor(fcen, df, nf, mode_reg)
    # the LDOS is cheap, so it is sampled finely enough (401 points) to
    # resolve the loaded-cavity Lorentzian
    ldos = mp.Ldos(fcen, df, 401)

    # ring-down for Q
    h = mp.Harminv(mp.Ey, src_pt, fcen, df)

    # effective mode volume at the dipole position from the ringing
    # field, V = int eps|E|^2 dV / (eps(r0) |E_y(r0)|^2), both averaged
    # over four optical periods starting 60 time units after the source
    # has switched off (transients gone, only the cavity mode rings).
    # Unlike the LDOS, this does not require the fields to decay fully,
    # so it gives the Purcell factor also for high-Q loadings.
    inner = mp.Volume(center=mp.Vector3(0.5 * (xin_L + xin_R), 0, zc) + shift,
                      size=mp.Vector3(xin_R - xin_L, 2 * yin, zt - zb))
    vs = dict(U=0.0, E2=0.0, n=0)

    def sample_volume(sim_):
        t_src_end = sim_.fields.last_source_time()
        t = sim_.meep_time()
        t1 = t_src_end + 60.0
        if t1 <= t < t1 + 4.0 / fcen:
            vs["U"] += sim_.electric_energy_in_box(box=inner)
            vs["E2"] += abs(sim_.get_field_point(mp.Ey, src_pt)) ** 2
            vs["n"] += 1

    # fixed ring-down time t_after after the source (a field-decay
    # stopping criterion is not met because weakly excited high-Q modes
    # elsewhere in the band keep ringing).  For N_out = 0 (Q_L ~ 900,
    # energy-decay time ~90 time units) t_after = 1000 lets the field
    # decay completely, so the LDOS spectrum is exact; for larger N_out a
    # shorter ring-down suffices for Q (Harminv), the mode volume and the
    # power fractions, which are ratios of the same ringing field
    sim.run(mp.dft_ldos(ldos=ldos), mp.after_sources(h), sample_volume,
            until_after_sources=t_after)
    eps0 = N_DIA ** 2
    V_um3 = 2.0 * vs["U"] / (eps0 * vs["E2"]) if vs["n"] else float("nan")
    freqs = np.array(mp.get_flux_freqs(f_right))
    P = {k: np.array(mp.get_fluxes(m)) for k, m in
         (("right", f_right), ("left", f_left), ("up", f_up), ("down", f_dn),
          ("yp", f_yp), ("ym", f_ym))}
    res = sim.get_eigenmode_coefficients(f_mode, [1], eig_parity=mp.ODD_Y,
                                         direction=mp.X)
    a_fwd = res.alpha[0, :, 0]
    P_te0 = np.abs(a_fwd) ** 2
    ld = np.array(sim.ldos_data)
    modes = [(m.freq, m.Q, abs(m.amp)) for m in h.modes]
    out = dict(n_out=n_out, resolution=resolution, freqs=freqs.tolist(),
               flux={k: v.tolist() for k, v in P.items()},
               P_te0=P_te0.tolist(), ldos=ld.tolist(), harminv=modes,
               V_um3=V_um3, V_samples=vs["n"], t_after=t_after,
               runtime_s=time.time() - t0, cell=[sx, sy, sz])
    return out


def bulk_ldos(resolution=40, fcen=1.0 / LAM, df_rel=0.05, nf=101, dpml=0.5):
    """LDOS of the same dipole in homogeneous diamond (Purcell reference)."""
    df = df_rel * fcen
    s = 2.0 + 2 * dpml
    sim = mp.Simulation(cell_size=mp.Vector3(s, s, s), resolution=resolution,
                        default_material=mp.Medium(index=N_DIA),
                        boundary_layers=[mp.PML(dpml)],
                        sources=[mp.Source(mp.GaussianSource(fcen, fwidth=df),
                                           component=mp.Ey, center=mp.Vector3())],
                        symmetries=[mp.Mirror(mp.Y, phase=-1), mp.Mirror(mp.X),
                                    mp.Mirror(mp.Z)])
    ldos = mp.Ldos(fcen, df, nf)
    sim.run(mp.dft_ldos(ldos=ldos), until_after_sources=mp.stop_when_fields_decayed(
        10, mp.Ey, mp.Vector3(), 1e-6))
    return dict(freqs=list(np.linspace(fcen - df / 2, fcen + df / 2, nf)),
                ldos=list(np.array(sim.ldos_data)))


if __name__ == "__main__":
    tag = sys.argv[1]
    res = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    if tag == "bulk":
        out = bulk_ldos(res)
        (RES_DIR / f"cavity_fdtd_bulk_r{res}.json").write_text(json.dumps(out))
    else:
        n_out = int(tag)
        ta = float(sys.argv[3]) if len(sys.argv) > 3 else 450.0
        out = run(n_out, res, t_after=ta)
        (RES_DIR / f"cavity_fdtd_N{n_out}_r{res}.json").write_text(json.dumps(out))
    print("done", tag, res)
