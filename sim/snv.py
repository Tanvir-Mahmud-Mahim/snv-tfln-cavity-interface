"""Strained SnV- spin Hamiltonian: eigenstructure, cyclicity, Rabi rate.

Implements the effective ground/excited-manifold Hamiltonian of the
negatively charged tin-vacancy centre in diamond exactly as written in
Rosenthal et al., Phys. Rev. X 13, 031022 (2023), Appendix B
(arXiv:2306.13199), Eqs. (B1)-(B6) and (B10)-(B17):

    H^{g,e} = H_SO + H_JT + H_Z + H_L                          (B5)

in the basis {|ex up>, |ex dn>, |ey up>, |ey dn>} = orbital (x/y) x spin:

    H_SO = -(hbar lambda/2) [[0, i], [-i, 0]] (x) sigma_z       (B1)
    H_JT =  hbar [[Ux, Uy], [Uy, -Ux]] (x) 1                    (B2)
    H_Z  = (hbar gamma/2) 1 (x) [[(1+2d)Bpar, Bperp],
                                 [Bperp*, -(1+2d)Bpar]]         (B3)
    H_L  = (hbar gamma f/2) [[0, i Bpar], [-i Bpar, 0]] (x) 1   (B4)

with Bpar = Bz, Bperp = Bx + i By in the spin frame (z along the SnV
high-symmetry axis <111>), gamma/2pi = 28 GHz/T, quenching factor f and
g-factor anisotropy delta given by f = gL*p, delta = gL*deltap (B16-B17).

Optical dipole operators between the manifolds (unpolarised light),
Eq. (B12):  px = q sigma_z, py = -q sigma_x, pz = q 1 (orbital space),
identity in spin.  Branching ratio (cyclicity of A1 vs A2, and likewise
C1 vs C2):  Lambda = P_1/P_2 with P_i = |<e_low| p.E |g_i>|^2 summed
over polarisation (B13)-(B15).

Parameter values: Table I of the same paper (device values):
  Delta_g/2pi = 902.98 GHz, Ups_g/2pi = 177.67 GHz, lam_g/2pi = 830.15 GHz
  Delta_e/2pi = 3000 GHz,   Ups_e/2pi = 134.00 GHz, lam_e/2pi = 2988.0 GHz
  gL_g = 0.363, gL_e = 0.581 -> f_g = 0.171, f_e = 0.073,
  delta_g = 0.0152, delta_e = 0.1761;  theta = 125.3 deg (dipole polar
  angle, fixed), phi = 37.33 deg (fit); ab initio reduction parameters
  deltap_g = 0.042, deltap_e = 0.303, p_g = 0.471, p_e = 0.125
  (Thiering & Gali, Phys. Rev. X 8, 021063 (2018)).

Validation targets (measured):
  * Rosenthal et al. 2023: Delta_g/2pi = 903 GHz (PL);
    qubit frequency vs field-angle sweep; Rabi ~ MHz-scale at b ~ mT.
  * Rosenthal et al., arXiv:2403.13110: at |B| = 180 mT cyclicity
    Lambda = 2244 +/- 108 at zeta = 147 deg and 8.6 +/- 0.4 at
    zeta = 53 deg; qubit frequency 3.677 GHz at |B| = 125 mT;
    Rabi frequency 6.25 MHz (80 ns pi pulse) at b ~ 0.6 mT.

All frequencies are angular frequencies over 2pi in GHz (i.e. H/(h)) so
eigenvalues read directly in GHz.
"""

from __future__ import annotations

import numpy as np
import qutip as qt

# ---------------------------------------------------------------- constants
GAMMA = 28.0        # electron gyromagnetic ratio, GHz/T  (gamma/2pi)

# Table I (Rosenthal et al., PRX 13, 031022 (2023)) -- device values
PARAMS = dict(
    lam_g=830.15, ups_g=177.67,        # GHz
    lam_e=2988.0, ups_e=134.00,        # GHz
    gL_g=0.363, gL_e=0.581,
    deltap_g=0.042, deltap_e=0.303, p_g=0.471, p_e=0.125,
    theta=np.deg2rad(125.3), phi=np.deg2rad(37.33),
)
PARAMS["f_g"] = PARAMS["gL_g"] * PARAMS["p_g"]          # 0.171
PARAMS["f_e"] = PARAMS["gL_e"] * PARAMS["p_e"]          # 0.073
PARAMS["delta_g"] = PARAMS["gL_g"] * PARAMS["deltap_g"]  # 0.0152
PARAMS["delta_e"] = PARAMS["gL_e"] * PARAMS["deltap_e"]  # 0.1761

SI = qt.qeye(2)
SX, SY, SZ = qt.sigmax(), qt.sigmay(), qt.sigmaz()


def h_manifold(lam, ups_x, ups_y, f, delta, B_spin):
    """4x4 manifold Hamiltonian (GHz), Eqs. B1-B5.

    B_spin = (Bx, By, Bz) in tesla, in the spin frame (z along <111>).
    """
    Bx, By, Bz = B_spin
    h_so = -0.5 * lam * qt.tensor(qt.Qobj([[0, 1j], [-1j, 0]]), SZ)
    h_jt = qt.tensor(qt.Qobj([[ups_x, ups_y], [ups_y, -ups_x]]), SI)
    zee = qt.Qobj([[(1 + 2 * delta) * Bz, Bx - 1j * By],
                   [Bx + 1j * By, -(1 + 2 * delta) * Bz]])
    h_z = 0.5 * GAMMA * qt.tensor(SI, zee)
    h_l = 0.5 * GAMMA * f * Bz * qt.tensor(
        qt.Qobj([[0, 1j], [-1j, 0]]), SI)
    return h_so + h_jt + h_z + h_l


def lab_to_spin(B_lab, theta, phi):
    """Rotate a lab-frame vector into the spin frame whose z-axis points
    along the dipole direction mu(theta, phi) (Fig. 9 of Rosenthal 2023).

    The spin frame is built by rotating the lab frame with R_z(phi) then
    R_y(theta); vectors transform with the inverse rotation.
    """
    ct, st = np.cos(theta), np.sin(theta)
    cp, sp = np.cos(phi), np.sin(phi)
    # rows are the spin-frame basis vectors expressed in the lab frame
    ex = np.array([ct * cp, ct * sp, -st])
    ey = np.array([-sp, cp, 0.0])
    ez = np.array([st * cp, st * sp, ct])       # dipole direction
    R = np.vstack([ex, ey, ez])
    return R @ np.asarray(B_lab)


def b_lab(zeta_deg, Bmag):
    """Magnet field in the lab frame: swept on the X/Z circle (Fig. 9),
    B = |B| (cos zeta, 0, sin zeta)."""
    z = np.deg2rad(zeta_deg)
    return Bmag * np.array([np.cos(z), 0.0, np.sin(z)])


def solve(B_lab_vec, prm=PARAMS, ups_scale=1.0):
    """Eigenstates of ground and excited manifolds.

    Returns (evals_g, ekets_g, evals_e, ekets_e), energies in GHz sorted
    ascending.  ups_scale scales the strain of both manifolds (strain
    engineering axis of the paper).  The strain axis is taken along x
    (Ups_y = 0), which fixes the azimuthal origin of the spin frame.
    """
    Bs = lab_to_spin(B_lab_vec, prm["theta"], prm["phi"])
    hg = h_manifold(prm["lam_g"], ups_scale * prm["ups_g"], 0.0,
                    prm["f_g"], prm["delta_g"], Bs)
    he = h_manifold(prm["lam_e"], ups_scale * prm["ups_e"], 0.0,
                    prm["f_e"], prm["delta_e"], Bs)
    eg, kg = hg.eigenstates()
    ee, ke = he.eigenstates()
    return eg, kg, ee, ke


# dipole operators (orbital part, Eq. B12); identity in spin
PX = qt.tensor(qt.Qobj([[1, 0], [0, -1]]), SI)
PY = qt.tensor(qt.Qobj([[0, -1], [-1, 0]]), SI)
PZ = qt.tensor(qt.Qobj([[1, 0], [0, 1]]), SI)


def transition_strengths(ket_g, ket_e):
    """|<e|p_i|g>|^2 summed over the three dipole components
    (unpolarised light, Eq. B13/B14)."""
    return sum(abs(complex(ket_e.dag() * P * ket_g)) ** 2
               for P in (PX, PY, PZ))


def cyclicity(B_lab_vec, prm=PARAMS, ups_scale=1.0, manifold="lower"):
    """Cyclicity Lambda = P(spin-conserving)/P(spin-flipping) for the
    transitions between the two lowest ground states |1>,|2> and the
    relevant excited state (Eq. B15).

    manifold='lower' uses the lowest excited doublet (C/D transitions,
    states |C>,|D> at index 0,1), 'upper' the A/B doublet (index 2,3).
    The spin-conserving partner of |1> is whichever excited state has
    the larger dipole overlap.
    """
    eg, kg, ee, ke = solve(B_lab_vec, prm, ups_scale)
    i0 = 0 if manifold == "lower" else 2
    exc = ke[i0], ke[i0 + 1]
    # A1-like: |1> -> spin-conserving partner; A2-like: |2> -> same state
    s11 = transition_strengths(kg[0], exc[0])
    s12 = transition_strengths(kg[0], exc[1])
    if s11 >= s12:
        e_cyc = exc[0]
    else:
        e_cyc = exc[1]
    p_conserve = transition_strengths(kg[0], e_cyc)
    p_flip = transition_strengths(kg[1], e_cyc)
    return p_conserve / p_flip


def qubit_frequency(B_lab_vec, prm=PARAMS, ups_scale=1.0):
    eg, kg, ee, ke = solve(B_lab_vec, prm, ups_scale)
    return float(eg[1] - eg[0])


def rabi_rate(B_lab_vec, b_lab_vec, prm=PARAMS, ups_scale=1.0):
    """Microwave Rabi rate Omega_MW/2pi (GHz) = |<2|H_MW|1>| (Eq. B10-B11),
    with the drive field b given in the lab frame in tesla."""
    eg, kg, ee, ke = solve(B_lab_vec, prm, ups_scale)
    bs = lab_to_spin(b_lab_vec, prm["theta"], prm["phi"])
    bx, by, bz = bs
    zee = qt.Qobj([[bz, bx - 1j * by], [bx + 1j * by, -bz]])
    h_mw = 0.5 * GAMMA * qt.tensor(SI, zee)
    return abs(complex(kg[1].dag() * h_mw * kg[0]))


if __name__ == "__main__":
    p = PARAMS
    dg = np.sqrt(p["lam_g"] ** 2 + 4 * p["ups_g"] ** 2)
    de = np.sqrt(p["lam_e"] ** 2 + 4 * p["ups_e"] ** 2)
    print(f"Delta_g = {dg:.2f} GHz (target 902.98 +/- 0.73)")
    print(f"Delta_e = {de:.1f} GHz (target 3000)")
    print(f"f_g={p['f_g']:.3f} f_e={p['f_e']:.3f} "
          f"delta_g={p['delta_g']:.4f} delta_e={p['delta_e']:.4f}")

    # zero-field splitting from full diagonalisation
    eg, kg, ee, ke = solve([0, 0, 0])
    print(f"diag: Delta_g = {eg[2]-eg[0]:.2f} GHz, "
          f"Delta_e = {ee[2]-ee[0]:.1f} GHz")

    # qubit frequency at the readout-paper operating point
    for zeta, B, target in ((147.0, 0.125, 3.677),):
        fq = qubit_frequency(b_lab(zeta, B))
        print(f"qubit freq at zeta={zeta} deg, B={B*1e3:.0f} mT: "
              f"{fq:.3f} GHz (measured {target})")

    # cyclicity vs zeta at 180 mT
    for zeta, target in ((147.0, 2244.0), (53.0, 8.6)):
        lam = cyclicity(b_lab(zeta, 0.180))
        print(f"cyclicity at zeta={zeta:5.1f} deg, 180 mT: {lam:9.1f} "
              f"(measured {target})")

    # Rabi rate at 0.6 mT drive (perpendicular to dipole), B=125 mT
    Bv = b_lab(147.0, 0.125)
    for ang in (0, 45, 90):
        a = np.deg2rad(ang)
        bvec = 6e-4 * np.array([np.sin(a), 0, np.cos(a)])
        om = rabi_rate(Bv, bvec) * 1e3
        print(f"Rabi at |b|=0.6 mT, drive angle {ang:3d} deg: "
              f"{om:6.2f} MHz (measured 6.25)")
