"""FEM waveguide cross-section studies for the diamond-on-TFLN interface.

All mode solving uses femwell (vector FEM Maxwell mode solver on scikit-fem)
with gmsh meshes.  Geometry follows the heterogeneous diamond / thin-film
lithium-niobate platform of Riedel et al., "Efficient photonic integration
of diamond color centers and thin-film lithium niobate," (2023),
arXiv:2306.15207 / doi:10.1364/OPTICAQ.1.000103: diamond nanobeams of
~200 nm thickness with tapered ends (350 nm -> 50 nm over ~10 um in their
SiV design) placed on a LN-on-insulator chip with a 190-nm LN film on SiO2.
Here the target wavelength is the SnV zero-phonon line at 619 nm instead of
the SiV line at 737 nm, so widths are re-optimised for single-mode
operation at 619 nm.

Studies
-------
A.  Diamond nanobeam (t = 200 nm) on SiO2: n_eff versus width -> choice of
    the single-mode cavity-waveguide width.
B.  TFLN ridge (t = 190 nm) on SiO2: n_eff versus width -> choice of the
    LN bus-waveguide width.
C.  Hybrid cross-section (diamond beam on top of the LN ridge): supermode
    effective indices versus diamond width -> adiabatic-transfer
    anticrossing, and eigenmode-expansion (EME) transfer efficiency of the
    linear taper versus taper length.

Units: um throughout the geometry; wavelength in um.
"""

from __future__ import annotations

import json
import pathlib
from collections import OrderedDict

import numpy as np
from shapely.geometry import box

from skfem import Basis, ElementTriP0
from femwell.mesh import mesh_from_OrderedDict
from femwell.maxwell.waveguide import compute_modes, calculate_overlap

import materials

RESULTS = pathlib.Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

LAM = materials.LAMBDA_SNV_UM      # 0.619 um
N_DIA = float(materials.diamond.n(LAM))
N_LN_O = float(materials.linbo3_o.n(LAM))   # in-plane (x-cut film, TE sees n_o)
N_LN_E = float(materials.linbo3_e.n(LAM))
N_SIO2 = float(materials.silica.n(LAM))

T_DIA = 0.200    # diamond nanobeam thickness (Riedel et al.: ~200 nm)
T_LN = 0.190     # LN film thickness (Riedel et al.: 190 nm)


def _solve(shapes: OrderedDict, eps_map: dict, resolutions: dict,
           num_modes: int = 2, order: int = 1):
    mesh = mesh_from_OrderedDict(
        shapes, resolutions, default_resolution_max=0.05, verbose=False
    )
    from skfem.io.meshio import from_meshio
    m = from_meshio(mesh)
    basis0 = Basis(m, ElementTriP0())
    eps = basis0.zeros(dtype=complex)
    for name, val in eps_map.items():
        eps[basis0.get_dofs(elements=name)] = val
    return compute_modes(basis0, eps, wavelength=LAM, num_modes=num_modes,
                         order=order)


def diamond_on_silica(w_dia: float, num_modes: int = 3):
    """Diamond nanobeam (w x 200 nm) resting on SiO2, air above."""
    pad, hclad = 0.9, 0.8
    core = box(-w_dia / 2, 0, w_dia / 2, T_DIA)
    sub = box(-w_dia / 2 - pad, -hclad, w_dia / 2 + pad, 0)
    air = box(-w_dia / 2 - pad, 0, w_dia / 2 + pad, T_DIA + hclad)
    shapes = OrderedDict(core=core, sub=sub, air=air)
    res = {"core": {"resolution": 0.012, "distance": 0.4}}
    eps = {"core": N_DIA ** 2, "sub": N_SIO2 ** 2, "air": 1.0}
    return _solve(shapes, eps, res, num_modes=num_modes)


def ln_ridge(w_ln: float, num_modes: int = 3):
    """Fully-etched TFLN ridge (w x 190 nm) on SiO2, air above.

    x-cut film: TE (in-plane E) sees the ordinary index at propagation
    along y; we use n_o for the film in this scalar-isotropic treatment
    and quote the o/e values in the manuscript.
    """
    pad, hclad = 1.0, 0.8
    core = box(-w_ln / 2, 0, w_ln / 2, T_LN)
    sub = box(-w_ln / 2 - pad, -hclad, w_ln / 2 + pad, 0)
    air = box(-w_ln / 2 - pad, 0, w_ln / 2 + pad, T_LN + hclad)
    shapes = OrderedDict(core=core, sub=sub, air=air)
    res = {"core": {"resolution": 0.012, "distance": 0.4}}
    eps = {"core": N_LN_O ** 2, "sub": N_SIO2 ** 2, "air": 1.0}
    return _solve(shapes, eps, res, num_modes=num_modes)


def hybrid(w_dia: float, w_ln: float, num_modes: int = 4, order: int = 1,
           res_core: float = 0.012):
    """Diamond beam (w_dia x 200 nm) sitting on the LN ridge
    (w_ln x 190 nm) on SiO2."""
    pad, hclad = 1.0, 0.8
    half = max(w_dia, w_ln) / 2 + pad
    dia = box(-w_dia / 2, T_LN, w_dia / 2, T_LN + T_DIA)
    ln = box(-w_ln / 2, 0, w_ln / 2, T_LN)
    sub = box(-half, -hclad, half, 0)
    air = box(-half, 0, half, T_LN + T_DIA + hclad)
    shapes = OrderedDict(dia=dia, ln=ln, sub=sub, air=air)
    res = {"dia": {"resolution": res_core, "distance": 0.4},
           "ln": {"resolution": res_core, "distance": 0.4}}
    eps = {"dia": N_DIA ** 2, "ln": N_LN_O ** 2, "sub": N_SIO2 ** 2,
           "air": 1.0}
    return _solve(shapes, eps, res, num_modes=num_modes, order=order)


def te_modes(modes, min_te=0.7):
    return [m for m in modes if m.te_fraction > min_te]


# ----------------------------------------------------------------- studies

def study_A(widths=None):
    widths = widths if widths is not None else np.arange(0.16, 0.42, 0.02)
    out = {"w": [], "neff": [], "te": []}
    for w in widths:
        modes = diamond_on_silica(float(w))
        ne = [float(np.real(m.n_eff)) for m in modes]
        te = [float(m.te_fraction) for m in modes]
        out["w"].append(float(w)); out["neff"].append(ne); out["te"].append(te)
        print(f"  w={w*1e3:4.0f} nm  n_eff={['%.4f' % x for x in ne]}"
              f"  TE={['%.2f' % x for x in te]}", flush=True)
    (RESULTS / "wg_diamond.json").write_text(json.dumps(out, indent=1))
    return out


def study_B(widths=None):
    widths = widths if widths is not None else np.arange(0.25, 0.75, 0.05)
    out = {"w": [], "neff": [], "te": []}
    for w in widths:
        modes = ln_ridge(float(w))
        ne = [float(np.real(m.n_eff)) for m in modes]
        te = [float(m.te_fraction) for m in modes]
        out["w"].append(float(w)); out["neff"].append(ne); out["te"].append(te)
        print(f"  w={w*1e3:4.0f} nm  n_eff={['%.4f' % x for x in ne]}"
              f"  TE={['%.2f' % x for x in te]}", flush=True)
    (RESULTS / "wg_ln.json").write_text(json.dumps(out, indent=1))
    return out


def study_C_anticrossing(w_ln: float, widths=None):
    widths = widths if widths is not None else np.arange(0.05, 0.36, 0.02)
    out = {"w_ln": w_ln, "w": [], "neff": [], "te": []}
    for w in widths:
        modes = hybrid(float(w), w_ln)
        ne = [float(np.real(m.n_eff)) for m in modes]
        te = [float(m.te_fraction) for m in modes]
        out["w"].append(float(w)); out["neff"].append(ne); out["te"].append(te)
        print(f"  w_dia={w*1e3:4.0f} nm  n_eff={['%.4f' % x for x in ne]}",
              flush=True)
    (RESULTS / "wg_hybrid.json").write_text(json.dumps(out, indent=1))
    return out


# ------------------------------------------------------- EME taper transfer

def _taper_shared_mesh(ws, w_ln: float):
    """Build ONE mesh whose subdomains are nested diamond-width rings, so
    every taper cross-section can be represented on the same mesh (only the
    per-subdomain permittivity changes).  Required because femwell's
    calculate_overlap needs both modes on the same mesh."""
    pad, hclad = 1.0, 0.8
    half = max(ws[0], w_ln) / 2 + pad
    shapes = OrderedDict()
    # innermost diamond core (always diamond)
    shapes["core"] = box(-ws[-1] / 2, T_LN, ws[-1] / 2, T_LN + T_DIA)
    for i in range(len(ws) - 1):
        wo, wi = ws[i], ws[i + 1]
        shapes[f"ringL{i}"] = box(-wo / 2, T_LN, -wi / 2, T_LN + T_DIA)
        shapes[f"ringR{i}"] = box(wi / 2, T_LN, wo / 2, T_LN + T_DIA)
    shapes["ln"] = box(-w_ln / 2, 0, w_ln / 2, T_LN)
    shapes["sub"] = box(-half, -hclad, half, 0)
    shapes["air"] = box(-half, 0, half, T_LN + T_DIA + hclad)
    res = {"core": {"resolution": 0.015, "distance": 0.3},
           "ln": {"resolution": 0.015, "distance": 0.3}}
    for i in range(len(ws) - 1):
        res[f"ringL{i}"] = {"resolution": 0.015, "distance": 0.1}
        res[f"ringR{i}"] = {"resolution": 0.015, "distance": 0.1}
    mesh = mesh_from_OrderedDict(shapes, res, default_resolution_max=0.05,
                                 verbose=False)
    from skfem.io.meshio import from_meshio
    m = from_meshio(mesh)
    basis0 = Basis(m, ElementTriP0())
    return m, basis0


def _taper_eps(basis0, k: int, n_ring: int):
    """Permittivity vector for taper segment k (diamond down to ring k)."""
    eps = basis0.zeros(dtype=complex) + 1.0
    eps[basis0.get_dofs(elements="air")] = 1.0
    eps[basis0.get_dofs(elements="sub")] = N_SIO2 ** 2
    eps[basis0.get_dofs(elements="ln")] = N_LN_O ** 2
    eps[basis0.get_dofs(elements="core")] = N_DIA ** 2
    for i in range(n_ring):
        val = N_DIA ** 2 if i >= k else 1.0
        eps[basis0.get_dofs(elements=f"ringL{i}")] = val
        eps[basis0.get_dofs(elements=f"ringR{i}")] = val
    return eps


def eme_taper(w_ln: float, w0: float = 0.35, w1: float = 0.05,
              lengths=None, n_seg: int = 30, num_modes: int = 3):
    """Eigenmode-expansion transfer through a linear diamond taper
    (width w0 -> w1) on the LN ridge.

    The taper is discretised into n_seg segments; in each segment the local
    supermodes are computed, propagated with exp(i beta L_seg) and projected
    onto the modes of the next segment via the unconjugated field overlap
    (femwell.calculate_overlap).  Radiation loss beyond the retained mode
    set is treated as loss (non-unitary projection).  Launch: fundamental
    (diamond-like) supermode at w0.  Readout: power in the fundamental
    (LN-like) supermode at w1.  All segments share one mesh (nested-ring
    subdomains) so that modal overlaps are computed exactly.
    """
    lengths = lengths if lengths is not None else [2, 4, 6, 8, 10, 14, 18]
    ws = np.linspace(w0, w1, n_seg)
    print(f"  building shared mesh ({n_seg} rings) ...", flush=True)
    m, basis0 = _taper_shared_mesh(ws, w_ln)
    print(f"  mesh: {m.p.shape[1]} vertices; solving {n_seg} cross-sections",
          flush=True)
    packs = []
    for i in range(n_seg):
        eps = _taper_eps(basis0, i, n_seg - 1)
        modes = compute_modes(basis0, eps, wavelength=LAM,
                              num_modes=num_modes, order=1)
        packs.append(compute_pack(modes))
        print(f"    [{i+1}/{n_seg}] w={ws[i]*1e3:5.1f} nm  "
              f"n_eff={['%.4f' % p['neff'] for p in packs[-1]]}", flush=True)

    # unconjugated overlap matrices between adjacent segments
    overlaps = []
    for a, b in zip(packs[:-1], packs[1:]):
        M = np.zeros((len(b), len(a)), dtype=complex)
        for j, mb in enumerate(b):
            for i_, ma in enumerate(a):
                M[j, i_] = calculate_overlap(
                    ma["basis"], ma["E"], ma["H"], mb["basis"], mb["E"], mb["H"]
                )
        overlaps.append(M)

    res = {"w_ln": w_ln, "w0": w0, "w1": w1, "n_seg": n_seg,
           "lengths": list(map(float, lengths)), "T": []}
    for L in lengths:
        Lseg = L / n_seg
        amp = np.zeros(len(packs[0]), dtype=complex)
        amp[0] = 1.0          # launch fundamental supermode (diamond-like)
        for k, M in enumerate(overlaps):
            beta = np.array([2 * np.pi * p["neff"] / LAM for p in packs[k]])
            amp = M @ (amp * np.exp(1j * beta * Lseg))
        beta = np.array([2 * np.pi * p["neff"] / LAM for p in packs[-1]])
        amp = amp * np.exp(1j * beta * (L / n_seg))
        T = float(np.abs(amp[0]) ** 2)
        res["T"].append(T)
        print(f"  L={L:5.1f} um  T(fundamental)={T:.4f}", flush=True)
    (RESULTS / "taper_eme.json").write_text(json.dumps(res, indent=1))
    return res


def compute_pack(modes):
    """Extract basis/fields, normalised to unit self-overlap (power),
    sorted by n_eff."""
    packs = []
    for m in modes:
        P = calculate_overlap(m.basis, m.E, m.H, m.basis, m.E, m.H)
        E = m.E / np.sqrt(P)
        H = m.H / np.sqrt(P)
        packs.append({"neff": float(np.real(m.n_eff)),
                      "basis": m.basis, "E": E, "H": H})
    packs.sort(key=lambda p: -p["neff"])
    return packs


if __name__ == "__main__":
    print(f"n(diamond)={N_DIA:.4f} n(LN,o)={N_LN_O:.4f} "
          f"n(LN,e)={N_LN_E:.4f} n(SiO2)={N_SIO2:.4f} at {LAM} um")
    print("Study A: diamond nanobeam on SiO2")
    study_A()
    print("Study B: TFLN ridge on SiO2")
    study_B()
