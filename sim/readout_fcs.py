"""Hamiltonian-consistent single-shot readout model (revision 1).

Replaces the two-level counting model of ``readout.py`` for all
predictions of the revised manuscript.  Everything that the readout
depends on is computed from the strained SnV- Hamiltonian of ``snv.py``
at the actual magnetic field (magnitude and direction) and strain:

* the four optical transitions of the C line between the two ground
  spin states {|dn>, |up>} (lower orbital branch) and the two excited
  spin states {|A>, |B>} (lower excited orbital branch), labelled as in
  Rosenthal et al., Phys. Rev. X 14, 041008 (2024): spin-conserving A1
  (|dn> <-> |A>, driven by the readout laser) and B2 (|up> <-> |B>),
  spin-flipping A2 (|up> <-> |A>) and B1 (|dn> <-> |B>);
* their frequencies (qubit frequency w_q, excited splitting w_e, and
  hence the A1-B2 detuning w_e - w_q and the other detunings);
* their polarisation-summed dipole strengths S_ge (Eq. B13 of
  Rosenthal et al., PRX 13, 031022 (2023)), which set both the relative
  Rabi frequencies of the readout laser on each transition and the
  decay branching ratios (hence the cyclicities Lambda_A = S_A1/S_A2 and
  Lambda_B = S_B2/S_B1).

The emitter is coupled to a broad cavity (kappa >> w_q, w_e), so all
four C-line transitions are Purcell-enhanced by the same factor; the D
line (detuned by the ground-state splitting Delta_g) is enhanced
according to the cavity Lorentzian.  Optional extra pure dephasing
(homogeneous broadening beyond the Purcell-enhanced lifetime limit),
quasi-static spectral diffusion (Gaussian common-mode shift of all
transitions, fixed during one readout window), and laser background are
included.

Counting statistics are computed EXACTLY from the four-level Lindblad
master equation by full counting statistics: the detected-photon
generating function G(z) = Tr[exp(L_z tau) rho_0], with every decay
channel's jump term multiplied by (1 - eta + eta z), is evaluated on the
unit circle and inverted by FFT.  This includes repeated spin flips,
dark-state excitation through all transitions, and all coherent
saturation effects; no rate-equation or adiabatic approximation is made.

Units: frequencies in GHz (cycles/ns), rates in 1/ns, times in ns.
"""

from __future__ import annotations

import numpy as np

import snv

TWO_PI = 2.0 * np.pi

# ----------------------------------------------------- emitter constants
TAU_0 = 4.5            # bulk excited-state lifetime, ns [Rosenthal 2024]
GAMMA_0 = 1.0 / TAU_0  # 1/ns
ETA_Q = 0.80           # quantum efficiency [Rosenthal 2024]
ETA_DW = 0.57          # Debye-Waller factor [Goerlitz 2020]
ETA_BR = 0.75          # C/D branching ratio [Lee 2026]
XI_POL = 2.0 / 3.0     # |C . e_cav|^2 for the 45-deg beam [Lee 2026]
N_DIA = 2.4137
C0 = 299792458.0
LAM_C_NM = 618.9       # C line wavelength used for the cavity frequency


# ------------------------------------------------------- level structure

def level_structure(B_lab_vec, ups_scale=1.0):
    """Frequencies and polarisation-summed dipole strengths of the four
    C-line transitions at field B (lab frame, tesla) and strain scale.

    Returns dict with
      S   2x2 array S[g, e], g in (dn, up), e in (A, B)
      nu  2x2 array of transition frequencies relative to A1 (GHz)
      wq  ground (qubit) splitting, we excited splitting (GHz)
      LamA, LamB cyclicities of |A> and |B>
      delta_g  ground-state orbital splitting (GHz), for the D line
    The bright state |dn> is the lower ground spin state; |A> is the
    excited state it is optically cycled to (the larger strength)."""
    eg, kg, ee, ke = snv.solve(B_lab_vec, ups_scale=ups_scale)
    ts = snv.transition_strengths
    S = np.array([[ts(kg[g], ke[e]) for e in (0, 1)] for g in (0, 1)])
    a = int(np.argmax(S[0]))
    b = 1 - a
    S = S[:, [a, b]]
    Eg = np.array([eg[0], eg[1]]) - eg[0]
    Ee = np.array([ee[a], ee[b]])
    nu = Ee[None, :] - Eg[:, None]
    nu = nu - nu[0, 0]
    return dict(S=S, nu=nu, wq=float(Eg[1]), we=float(Ee[1] - Ee[0]),
                LamA=float(S[0, 0] / S[1, 0]), LamB=float(S[1, 1] / S[0, 1]),
                delta_g=float(eg[2] - eg[0]))


# -------------------------------------------------------- cavity budget

def purcell_max(Q_L, V_rel):
    return 3.0 / (4.0 * np.pi ** 2) * Q_L / V_rel


def cavity_budget(Q_L, V_rel=0.68, xi_pos=0.5, delta_cd_ghz=903.0,
                  eta_wg=None, Q_i=1.39e6, T_taper=None, eta_chip=0.5,
                  eta_det=0.77, xi_pol_D=XI_POL, f_cav_ghz=None):
    """Transition-resolved Purcell budget of the C-line-resonant cavity.

    The cavity-mode emission rate of a transition with Purcell factor F
    is (F - 1) times its free-space rate (the free-space channel
    persists), so the fraction of all decays emitted into the cavity mode
    is beta = (zeta - 1)/zeta with
      zeta = 1 + eta_q eta_DW [eta_BR (F_C - 1) + (1 - eta_BR)(F_D - 1)].
    The D line is detuned from the cavity by the ground-state splitting
    Delta_g and enhanced according to the cavity Lorentzian."""
    f_cav = f_cav_ghz if f_cav_ghz is not None else C0 / (LAM_C_NM * 1e-9) / 1e9
    kappa = f_cav / Q_L                        # FWHM, GHz
    F_max = purcell_max(Q_L, V_rel)
    F_C = 1.0 + (F_max - 1.0) * XI_POL * xi_pos
    L_D = 1.0 / (1.0 + (2.0 * delta_cd_ghz / kappa) ** 2)
    F_D = 1.0 + (F_max - 1.0) * xi_pol_D * xi_pos * L_D
    zeta = 1.0 + ETA_Q * ETA_DW * (ETA_BR * (F_C - 1.0)
                                   + (1.0 - ETA_BR) * (F_D - 1.0))
    beta = (zeta - 1.0) / zeta
    if T_taper is None:
        raise ValueError("pass the EME taper transfer T_taper")
    if eta_wg is None:
        # closed-cavity estimate (only if no FDTD result is available)
        eta_wg = 1.0 - Q_L / Q_i
    eta = beta * eta_wg * T_taper * eta_chip * eta_det
    return dict(Q_L=Q_L, kappa_GHz=kappa, F_max=F_max, F_C=F_C, F_D=F_D,
                zeta=zeta, beta=beta, eta_wg=eta_wg, eta=eta,
                gamma_GHz=zeta * GAMMA_0 / TWO_PI,
                tau_ns=TAU_0 / zeta)


# ---------------------------------------------------- Lindblad machinery

def _superop_terms(ls, gamma_tot, s, gamma_phi_extra, delta_L):
    """Build L0 (everything except the jump terms) and J (sum of jump
    terms) for the 4-level system in the laser rotating frame.

    States: 0=|dn>, 1=|up>, 2=|A>, 3=|B>.  gamma_tot: total excited-state
    decay rate (1/ns).  s: saturation parameter of A1.
    gamma_phi_extra: extra homogeneous FWHM (GHz) beyond the lifetime
    limit, implemented as pure dephasing.  delta_L: laser detuning from
    A1 (GHz)."""
    S, nu = ls["S"], ls["nu"]
    d = 4
    gstar = TWO_PI * gamma_phi_extra / 2.0     # coherence dephasing, 1/ns
    G2 = gamma_tot / 2.0 + gstar
    Om0 = np.sqrt(s * gamma_tot * G2)           # A1 Rabi frequency (rad/ns)
    H = np.zeros((d, d), dtype=complex)
    H[1, 1] = TWO_PI * ls["wq"]
    H[2, 2] = -TWO_PI * delta_L
    H[3, 3] = TWO_PI * (nu[0, 1] - delta_L)
    for g in (0, 1):
        for e in (0, 1):
            Om = Om0 * np.sqrt(S[g, e] / S[0, 0])
            H[2 + e, g] = H[g, 2 + e] = Om / 2.0
    I = np.eye(d)
    L0 = -1j * (np.kron(I, H) - np.kron(H.T, I))
    J = np.zeros((d * d, d * d), dtype=complex)
    for e in (0, 1):
        Stot = S[0, e] + S[1, e]
        for g in (0, 1):
            c = np.zeros((d, d))
            c[g, 2 + e] = np.sqrt(gamma_tot * S[g, e] / Stot)
            cdc = c.T @ c
            J += np.kron(c.conj(), c)
            L0 += -0.5 * np.kron(I, cdc) - 0.5 * np.kron(cdc.T, I)
    if gstar > 0:
        P = np.zeros((d, d)); P[2, 2] = P[3, 3] = 1.0
        c = np.sqrt(2.0 * gstar) * P
        cdc = c.T @ c
        L0 += np.kron(c, c) - 0.5 * np.kron(I, cdc) - 0.5 * np.kron(cdc.T, I)
    return L0, J


def _vec(rho):
    return rho.reshape(-1, order="F")


def generating_functions(ls, gamma_tot, eta, s, taus, n_fft=1024,
                         gamma_phi_extra=0.0, delta_L=0.0, bg_rate=0.0):
    """Detected-photon pmfs P(k|dn) and P(k|up) for each window in taus.

    eta: detection probability per emitted photon (any transition).
    bg_rate: background count rate (1/ns), Poisson.
    Returns arrays Pb, Pd of shape (len(taus), n_fft)."""
    L0, J = _superop_terms(ls, gamma_tot, s, gamma_phi_extra, delta_L)
    zs = np.exp(2j * np.pi * np.arange(n_fft) / n_fft)
    # L_z = L0 + (1 - eta + eta z) J
    Lz = L0[None, :, :] + (1.0 - eta + eta * zs)[:, None, None] * J[None]
    w, V = np.linalg.eig(Lz)
    Vinv = np.linalg.inv(V)
    d = 4
    vecI = _vec(np.eye(d))
    r0 = np.zeros((d, d)); r0[0, 0] = 1.0
    r1 = np.zeros((d, d)); r1[1, 1] = 1.0
    left = np.einsum("i,zij->zj", vecI, V)               # (n_fft, 16)
    rb = np.einsum("zij,j->zi", Vinv, _vec(r0))
    rd = np.einsum("zij,j->zi", Vinv, _vec(r1))
    taus = np.atleast_1d(taus)
    E = np.exp(w[None, :, :] * taus[:, None, None])       # (nt, n_fft, 16)
    Gb = np.einsum("zj,tzj,zj->tz", left, E, rb)
    Gd = np.einsum("zj,tzj,zj->tz", left, E, rd)
    if bg_rate > 0:
        bgf = np.exp((zs[None, :] - 1.0) * bg_rate * taus[:, None])
        Gb = Gb * bgf
        Gd = Gd * bgf
    Pb = np.real(np.fft.fft(Gb, axis=1)) / n_fft
    Pd = np.real(np.fft.fft(Gd, axis=1)) / n_fft
    return np.clip(Pb, 0, None), np.clip(Pd, 0, None)


def best_threshold(Pb, Pd):
    """Optimal integer threshold: declare bright if K >= Nr.
    Returns (F, Nr) with F = 1 - (P(K>=Nr|dark) + P(K<Nr|bright))/2."""
    cb = np.cumsum(Pb)
    cd = np.cumsum(Pd)
    Nr = np.arange(1, len(Pb))
    err_b = cb[Nr - 1]               # P(K < Nr | bright)
    err_d = 1.0 - cd[Nr - 1]         # P(K >= Nr | dark)
    F = 1.0 - 0.5 * (err_b + err_d)
    i = int(np.argmax(F))
    return float(F[i]), int(Nr[i])


def gauss_hermite_fwhm(fwhm, n=9):
    """Nodes/weights for averaging over a Gaussian of given FWHM."""
    if fwhm <= 0:
        return np.array([0.0]), np.array([1.0])
    x, w = np.polynomial.hermite_e.hermegauss(n)
    sig = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return sig * x, w / w.sum()


def readout_fidelity(ls, gamma_tot, eta, s_list=(0.2, 0.5, 1.0, 2.0),
                     taus=None, n_fft=1024, gamma_phi_extra=0.0,
                     sd_fwhm=0.0, bg_dark=1e-7, sbr=np.inf, f0=0.99,
                     n_sd=9, return_all=False):
    """Single-shot readout fidelity optimised over drive s, window tau and
    integer threshold.  Spectral diffusion: the count distributions are
    averaged over a quasi-static Gaussian laser detuning (FWHM sd_fwhm,
    GHz) BEFORE the threshold is chosen (the threshold cannot adapt to
    the unknown shot-to-shot detuning).

    bg_dark: detector dark-count rate (1/ns); sbr: ratio of the detected
    bright photon rate at s = 1 to the detected laser-background rate at
    s = 1 (laser background scales linearly with s)."""
    if taus is None:
        taus = np.geomspace(2.0, 2.0e5, 40)      # 2 ns .. 200 us
    taus = np.asarray(taus, dtype=float)
    dls, wts = gauss_hermite_fwhm(sd_fwhm, n_sd)
    best = None
    for s in s_list:
        R1 = 0.5 * gamma_tot * 1.0 / 2.0           # A1 scattering rate at s=1
        bg = bg_dark + (0.0 if np.isinf(sbr) else eta * R1 * s / sbr)
        Pb = 0.0; Pd = 0.0
        for dl, wt in zip(dls, wts):
            b, d = generating_functions(ls, gamma_tot, eta, s, taus, n_fft,
                                        gamma_phi_extra, dl, bg)
            Pb = Pb + wt * b
            Pd = Pd + wt * d
        for it, tau in enumerate(taus):
            # skip windows whose photon numbers could alias in the FFT:
            # the emission rate never exceeds gamma_tot
            if (eta * gamma_tot + bg) * tau > 0.6 * n_fft:
                continue
            F, Nr = best_threshold(Pb[it], Pd[it])
            F = f0 * F + (1 - f0) * 0.5
            kk = np.arange(Pb.shape[1])
            nb = float((kk * Pb[it]).sum()); nd = float((kk * Pd[it]).sum())
            if best is None or F > best["F"]:
                best = dict(F=F, Nr=Nr, s=s, tau_ns=float(tau), nb=nb, nd=nd)
    return best


def fidelity_vs_tau(ls, gamma_tot, eta, taus,
                    s_list=(0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
                    n_fft=1024, gamma_phi_extra=0.0, sd_fwhm=0.0,
                    bg_dark=1e-7, sbr=np.inf, f0=0.99, n_sd=9):
    """Best fidelity (over drive s and threshold) at each fixed window.
    Returns arrays F[tau], s_opt[tau], Nr[tau], nb[tau], nd[tau]."""
    taus = np.asarray(taus, dtype=float)
    dls, wts = gauss_hermite_fwhm(sd_fwhm, n_sd)
    nt = len(taus)
    F = np.full(nt, np.nan); so = np.full(nt, np.nan)
    Nro = np.zeros(nt, int); nbo = np.full(nt, np.nan); ndo = np.full(nt, np.nan)
    kk = np.arange(n_fft)
    for s in s_list:
        R1 = 0.25 * gamma_tot
        bg = bg_dark + (0.0 if np.isinf(sbr) else eta * R1 * s / sbr)
        Pb = 0.0; Pd = 0.0
        for dl, wt in zip(dls, wts):
            b, d = generating_functions(ls, gamma_tot, eta, s, taus, n_fft,
                                        gamma_phi_extra, dl, bg)
            Pb = Pb + wt * b
            Pd = Pd + wt * d
        for it, tau in enumerate(taus):
            if (eta * gamma_tot + bg) * tau > 0.6 * n_fft:
                continue
            f, nr = best_threshold(Pb[it], Pd[it])
            f = f0 * f + (1 - f0) * 0.5
            if np.isnan(F[it]) or f > F[it]:
                F[it], so[it], Nro[it] = f, s, nr
                nbo[it] = float((kk * Pb[it]).sum())
                ndo[it] = float((kk * Pd[it]).sum())
    return F, so, Nro, nbo, ndo


def conditional_rates(ls, gamma_tot, s, gamma_phi_extra=0.0, delta_L=0.0):
    """Quasi-steady photon and spin-flip rates (1/ns) while the spin is in
    |dn> (bright) or |up> (dark), from the Lindblad steady state with the
    other ground state's population re-injected (valid for Lambda >> 1).
    Returns dict R_b, G_bd, R_d, G_db (photon emission rates and flip
    rates).  Used for interpretation only; fidelities use the exact
    counting statistics above."""
    S = ls["S"]
    out = {}
    for name, keep in (("b", 0), ("d", 1)):
        L0, J = _superop_terms(ls, gamma_tot, s, gamma_phi_extra, delta_L)
        d = 4
        I = np.eye(d)
        # rebuild with re-injection: every decay lands in `keep`
        Lr = L0.copy()
        for e in (0, 1):
            c = np.zeros((d, d))
            c[keep, 2 + e] = np.sqrt(gamma_tot)
            Lr += np.kron(c.conj(), c)
        # steady state: null vector of Lr
        w, V = np.linalg.eig(Lr)
        k = int(np.argmin(np.abs(w)))
        rho = V[:, k].reshape(d, d, order="F")
        rho = rho / np.trace(rho)
        pA, pB = np.real(rho[2, 2]), np.real(rho[3, 3])
        bA = S[:, 0] / S[:, 0].sum()
        bB = S[:, 1] / S[:, 1].sum()
        other = 1 - keep
        out["R_" + name] = gamma_tot * (pA + pB)
        out["G_" + name + ("d" if keep == 0 else "b")] = gamma_tot * (
            pA * bA[other] + pB * bB[other])
        out["pA_" + name] = pA
        out["pB_" + name] = pB
        out["flip_viaA_" + name] = gamma_tot * pA * bA[other]
        out["flip_viaB_" + name] = gamma_tot * pB * bB[other]
    return out
