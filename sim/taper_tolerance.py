"""Fabrication and placement tolerance of the diamond -> TFLN adiabatic
taper (revision 1, Reviewer 1 point 1).

Eigenmode expansion (EME) of the linear diamond taper (350 -> 50 nm,
200 nm thick) on the LN bus, generalised to non-ideal geometries:

  * lateral offset dx of the diamond taper from the LN ridge centre;
  * angular misalignment theta (the offset varies linearly along the
    taper, pivoting at the taper centre);
  * a vertical air gap g between the LN top surface and the diamond;
  * sloped LN sidewalls (angle psi from the horizontal, ridge TOP width
    fixed at w_top) and partial etch (ridge height h_r on a residual
    slab of thickness t_slab).

All taper cross-sections share ONE finite-element mesh whose diamond
layer is cut into narrow vertical strips (default 5 nm); each segment
assigns the diamond permittivity to the strips it covers (fractional
filling at the two edges, linear permittivity averaging).  Modal
overlaps between adjacent segments are therefore exact finite-element
integrals, exactly as in waveguide.eme_taper.  After the last taper
segment the field is projected onto the modes of the bare LN bus (the
diamond tip terminates), so the reported transfer T is the power
delivered to the fundamental TE mode of the LN waveguide alone.

Usage:  python taper_tolerance.py <tag>   (see CONFIGS below)
"""

from __future__ import annotations

import ctypes
import gc
import json
import pathlib
import sys
import time
from collections import OrderedDict

import numpy as np
from shapely.geometry import Polygon, box
from skfem import Basis, ElementTriP0
from skfem.io.meshio import from_meshio
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes

import materials
from waveguide import compute_pack, LAM, N_DIA, N_LN_O, N_SIO2

RESULTS = pathlib.Path(__file__).resolve().parent / "results" / "tolerance_eme"
RESULTS.mkdir(parents=True, exist_ok=True)

T_DIA = 0.200


def build_mesh(x_extent, t_film, h_ridge, psi_deg, w_top, gap, dx_strip,
               pad=1.0, hclad=0.8):
    """Shared mesh: strips of the diamond layer from -x_extent..x_extent,
    LN ridge (trapezoid) on optional slab, SiO2 below, air above."""
    t_slab = t_film - h_ridge
    tan_psi = np.tan(np.deg2rad(psi_deg)) if psi_deg < 89.999 else np.inf
    run = 0.0 if np.isinf(tan_psi) else h_ridge / tan_psi
    w_bot = w_top + 2 * run
    half = max(x_extent, w_bot / 2) + pad
    y0 = t_film + gap                        # diamond bottom
    edges = np.arange(-x_extent, x_extent + 1e-9, dx_strip)
    shapes = OrderedDict()
    for i in range(len(edges) - 1):
        shapes[f"s{i}"] = box(edges[i], y0, edges[i + 1], y0 + T_DIA)
    ridge = Polygon([(-w_bot / 2, t_slab), (w_bot / 2, t_slab),
                     (w_top / 2, t_film), (-w_top / 2, t_film)])
    shapes["ln"] = ridge
    if t_slab > 1e-6:
        shapes["slab"] = box(-half, 0.0, half, t_slab)
    shapes["sub"] = box(-half, -hclad, half, 0.0)
    shapes["air"] = box(-half, 0.0, half, y0 + T_DIA + hclad)
    res = {"ln": {"resolution": 0.015, "distance": 0.3}}
    for i in range(len(edges) - 1):
        res[f"s{i}"] = {"resolution": 0.015, "distance": 0.1}
    mesh = mesh_from_OrderedDict(shapes, res, default_resolution_max=0.05,
                                 verbose=False)
    m = from_meshio(mesh)
    basis0 = Basis(m, ElementTriP0())
    return m, basis0, edges


def eps_vector(basis0, edges, x_c, w, has_slab):
    eps = basis0.zeros(dtype=complex) + 1.0
    eps[basis0.get_dofs(elements="sub")] = N_SIO2 ** 2
    eps[basis0.get_dofs(elements="ln")] = N_LN_O ** 2
    if has_slab:
        eps[basis0.get_dofs(elements="slab")] = N_LN_O ** 2
    lo, hi = x_c - w / 2, x_c + w / 2
    for i in range(len(edges) - 1):
        a, b = edges[i], edges[i + 1]
        f = max(0.0, min(b, hi) - max(a, lo)) / (b - a)
        if f > 0:
            eps[basis0.get_dofs(elements=f"s{i}")] = f * N_DIA ** 2 + (1 - f)
    return eps


def run(tag, offset=0.0, tilt_deg=0.0, gap=0.0, psi_deg=90.0,
        t_film=0.190, h_ridge=0.190, w_top=0.600, lengths=(2, 3, 4, 6, 8, 10, 14),
        n_seg=30, num_modes=6, dx_strip=0.005, w0=0.35, w1=0.05):
    t0 = time.time()
    ws = np.linspace(w0, w1, n_seg)
    max_off = abs(offset) + (max(lengths) / 2) * np.tan(np.deg2rad(abs(tilt_deg)))
    x_extent = np.ceil((w0 / 2 + max_off + 0.01) / dx_strip) * dx_strip
    # a strip edge must not coincide with the ridge-top corner (gmsh then
    # fails on a degenerate segment when the gap is zero); shift by half a
    # strip in that case (only the 150-nm offset is affected)
    if np.any(np.abs(np.abs(np.arange(-x_extent, x_extent + 1e-9, dx_strip))
                     - w_top / 2) < 1e-6):
        x_extent += dx_strip / 2
    m, basis0, edges = build_mesh(x_extent, t_film, h_ridge, psi_deg, w_top,
                                  gap, dx_strip)
    has_slab = (t_film - h_ridge) > 1e-6
    print(f"[{tag}] mesh {m.p.shape[1]} vertices, {len(edges)-1} strips",
          flush=True)

    def _interp_modes(modes):
        """Normalised (unit power) transverse fields interpolated at the
        quadrature points; the basis objects are not kept (memory)."""
        packs = compute_pack(modes)
        bas = packs[0]["basis"]
        out = []
        for p in packs:
            E = np.asarray(bas.interpolate(p["E"])[0].value)
            H = np.asarray(bas.interpolate(p["H"])[0].value)
            out.append((p["neff"], E, H))
        return out, bas.dx

    # bare LN bus (tip removed) for the final projection
    eps_ln = eps_vector(basis0, edges, 0.0, 0.0, has_slab)
    ln_modes = compute_modes(basis0, eps_ln, wavelength=LAM,
                             num_modes=num_modes, order=1)
    ln_I, dxw = _interp_modes(ln_modes)
    te_ln = [i for i, mm in enumerate(sorted(ln_modes, key=lambda q: -np.real(q.n_eff)))
             if mm.te_fraction > 0.7]
    i_te0 = te_ln[0]
    del ln_modes

    def _cross(a, b):
        return a[0] * b[1] - a[1] * b[0]

    def _ovl(IA, IB):
        """Same integral as femwell.calculate_overlap on a shared mesh:
        0.5 * int (E_i^* x H_j + E_j x H_i^*) . z dA."""
        M = np.zeros((len(IB), len(IA)), dtype=complex)
        for j, (_, Eb, Hb) in enumerate(IB):
            for i_, (_, Ea, Ha) in enumerate(IA):
                M[j, i_] = 0.5 * np.sum((_cross(np.conj(Ea), Hb)
                                         + _cross(Eb, np.conj(Ha))) * dxw)
        return M

    def sweep(xcs):
        """Solve all segments sequentially, keeping only two in memory.
        Returns per-segment n_eff lists, adjacent overlaps, final
        projection onto the bare LN modes."""
        neffs, Ms = [], []
        prev = None
        for k in range(n_seg):
            eps = eps_vector(basis0, edges, xcs[k], ws[k], has_slab)
            modes = compute_modes(basis0, eps, wavelength=LAM,
                                  num_modes=num_modes, order=1)
            cur, _ = _interp_modes(modes)
            del modes
            # return freed solver memory to the system (long sweeps
            # otherwise grow past the memory of a small machine)
            gc.collect()
            try:
                ctypes.CDLL("libc.so.6").malloc_trim(0)
            except OSError:
                pass
            neffs.append([c[0] for c in cur])
            if prev is not None:
                Ms.append(_ovl(prev, cur))
            prev = cur
        Mf = _ovl(prev, ln_I)
        return neffs, Ms, Mf

    def propagate(neffs, Ms, Mf, L):
        Lseg = L / n_seg
        amp = np.zeros(len(neffs[0]), dtype=complex)
        amp[0] = 1.0
        for k, M in enumerate(Ms):
            beta = 2 * np.pi * np.array(neffs[k]) / LAM
            amp = M @ (amp * np.exp(1j * beta * Lseg))
        beta = 2 * np.pi * np.array(neffs[-1]) / LAM
        amp = amp * np.exp(1j * beta * Lseg)
        T_hyb = float(np.abs(amp[0]) ** 2)
        a_ln = Mf @ amp
        return T_hyb, float(np.abs(a_ln[i_te0]) ** 2), float(np.sum(np.abs(a_ln) ** 2))

    res = dict(tag=tag, offset_nm=offset * 1e3, tilt_deg=tilt_deg,
               gap_nm=gap * 1e3, psi_deg=psi_deg, t_film_nm=t_film * 1e3,
               h_ridge_nm=h_ridge * 1e3, w_top_nm=w_top * 1e3,
               n_seg=n_seg, num_modes=num_modes, dx_strip_nm=dx_strip * 1e3,
               lengths=[float(L) for L in lengths], T_hybrid=[], T_LN=[],
               P_retained=[], neff_start=None, neff_LN=[c[0] for c in ln_I])
    if tilt_deg == 0.0:
        xcs = np.full(n_seg, offset)
        neffs, Ms, Mf = sweep(xcs)
        res["neff_start"] = neffs[0]
        for L in lengths:
            th, tl, pr = propagate(neffs, Ms, Mf, L)
            res["T_hybrid"].append(th); res["T_LN"].append(tl)
            res["P_retained"].append(pr)
            print(f"[{tag}] L={L:5.1f} um  T_hyb={th:.4f}  T_LN(TE0)={tl:.4f}",
                  flush=True)
    else:
        for L in lengths:
            zs = np.linspace(0.0, L, n_seg)
            xcs = offset + (zs - L / 2) * np.tan(np.deg2rad(tilt_deg))
            neffs, Ms, Mf = sweep(xcs)
            th, tl, pr = propagate(neffs, Ms, Mf, L)
            res["T_hybrid"].append(th); res["T_LN"].append(tl)
            res["P_retained"].append(pr)
            print(f"[{tag}] L={L:5.1f} um  T_hyb={th:.4f}  T_LN(TE0)={tl:.4f}",
                  flush=True)
    res["runtime_s"] = time.time() - t0
    (RESULTS / f"{tag}.json").write_text(json.dumps(res, indent=1))
    return res


CONFIGS = {
    "nominal": dict(),
    "off25": dict(offset=0.025), "off50": dict(offset=0.050),
    "off75": dict(offset=0.075), "off100": dict(offset=0.100),
    "off150": dict(offset=0.150),
    "tilt0p5": dict(tilt_deg=0.5, lengths=(6,)), "tilt1": dict(tilt_deg=1.0, lengths=(6,)),
    "tilt2": dict(tilt_deg=2.0, lengths=(6,)), "tilt3": dict(tilt_deg=3.0, lengths=(6,)),
    "gap10": dict(gap=0.010), "gap20": dict(gap=0.020), "gap50": dict(gap=0.050),
    "sw80": dict(psi_deg=80.0), "sw70": dict(psi_deg=70.0), "sw60": dict(psi_deg=60.0),
    "rib60": dict(t_film=0.250, h_ridge=0.190),
    "rib60_sw70": dict(t_film=0.250, h_ridge=0.190, psi_deg=70.0),
    "combined": dict(offset=0.040, tilt_deg=0.5, gap=0.010, psi_deg=70.0,
                     lengths=(6,)),
    "nominal_m4": dict(num_modes=4, lengths=(6,)),
    "nominal_m8": dict(num_modes=8, lengths=(6,)),
    "nominal_seg45": dict(n_seg=45, lengths=(6,)),
    "nominal_dx2p5": dict(dx_strip=0.0025, lengths=(6,)),
}

if __name__ == "__main__":
    for tag in sys.argv[1:]:
        run(tag, **CONFIGS[tag])
