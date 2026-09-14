"""Cavity-enhanced single-shot readout model for the SnV- spin.

Combines
  (i)   the photonic design outputs (loaded Q, V, waveguide out-coupling
        efficiency from the GME cavity stage; taper transfer from the EME
        stage),
  (ii)  the SnV- emission budget (quantum efficiency, Debye-Waller
        factor, C/D branching ratio), and
  (iii) the two-level readout/photon-statistics model of Rosenthal et
        al., arXiv:2403.13110, Appendix A (Eqs. A6-A20),
to predict detected-photon numbers and single-shot readout fidelity of
the integrated diamond-on-TFLN spin-photon interface.

Emission budget (bulk emitter):
  gamma      total excited-state decay rate, 1/(4.5 ns)   [2403.13110]
  eta_q      radiative quantum efficiency, 0.8            [Table II ibid.]
  eta_DW     Debye-Waller (ZPL) fraction, 0.57            [Goerlitz et al.;
                                                           used by Lee et al.]
  eta_BR     C/D branching ratio, 0.75                    [Lee et al. 2511.05740]

Cavity coupling (transition-resolved, following Lee et al. Eqs. 1-10):
  F_C = F_max * xi_pol * xi_pos     Purcell factor of the C transition
  F_max = (3/(4 pi^2)) (lam/n)^3 Q_L/V
  xi_pol = |C_hat . e_cav|^2 = sin^2(theta_C) (1+sin 2 psi)/2 -> 2/3 at
           psi = 45 deg for <100>-oriented diamond (theta_C = 54.7 deg)
  xi_pos: spatial overlap of the emitter with the field antinode.

Lifetime-reduction (rate) factor with the cavity resonant on C1 and the
spin-flipping partner C2 detuned by the qubit frequency:
  zeta = 1 + eta_q eta_DW eta_BR (F_C - 1)                 (Lee Eq. 3 + QE)

Cyclicity enhancement: the cavity enhances only the spin-conserving C1
rate (C2 is detuned by omega_q >> kappa/2), so
  Lambda_cav = Lambda_0 * zeta / L(delta = omega_q)-correction,
computed below with the explicit Lorentzian factor.

Detected-photon budget per emission event:
  beta      = eta_q eta_DW eta_BR F_C / zeta_tot  (into the cavity mode)
  eta_wg    = kappa_wg / kappa_tot                (one-sided out-coupling)
  T_taper   diamond->TFLN adiabatic transfer      (EME stage; Riedel et
                                                   al. measured 0.92)
  eta_chip  on-chip routing + off-chip coupling
  eta_det   detector quantum efficiency

Readout statistics (Rosenthal Eqs. A8, A16-A20), threshold N_r = 1:
  nb - nd = eta (Lambda+1)(1 - exp(-Gamma_p tau)),  Gamma_p = R/(1+Lambda)
  F_r = 1/2 + (f0/2) (exp(-nd) - exp(-nb))
"""

from __future__ import annotations

import json
import pathlib

import numpy as np

RESULTS = pathlib.Path(__file__).resolve().parent / "results"
RESULTS.mkdir(exist_ok=True)

# ------------------------------------------------------ emitter constants
TAU_0 = 4.5e-9          # bulk lifetime, s [2403.13110: gamma = (4.5 ns)^-1]
GAMMA_0 = 1.0 / TAU_0   # total decay rate, 1/s
ETA_Q = 0.80            # quantum efficiency [2403.13110 Table II: 80-90%]
ETA_DW = 0.57           # Debye-Waller factor [Goerlitz 2020; Lee 2511.05740]
ETA_BR = 0.75           # C/D branching ratio [Lee 2511.05740]
XI_POL = np.sin(np.deg2rad(54.7)) ** 2   # = 2/3, beam at psi=45 deg
N_DIA = 2.4137          # diamond refractive index at 619 nm


def purcell_max(Q, V_rel):
    """Ideal Purcell factor, F = (3/4pi^2) (lam/n)^3 Q/V with V in units
    (lam/n)^3."""
    return 3.0 / (4 * np.pi ** 2) * Q / V_rel


def lorentzian_suppression(delta_hz, f0_hz, Q):
    """Purcell reduction for a transition detuned by delta from the
    cavity: L = 1 / (1 + (2 Q delta / f0)^2)."""
    return 1.0 / (1.0 + (2.0 * Q * delta_hz / f0_hz) ** 2)


class CavityInterface:
    """All derived quantities of the cavity-coupled SnV interface."""

    def __init__(self, Q_loaded, V_rel, eta_wg, xi_pos=0.8,
                 T_taper=0.95, eta_chip=0.8, eta_det=0.9,
                 Lambda_0=2244.0, omega_q_hz=3.677e9,
                 f_cav_hz=484.13e12):
        self.Q, self.V = Q_loaded, V_rel
        self.eta_wg = eta_wg
        self.xi_pos, self.T_taper = xi_pos, T_taper
        self.eta_chip, self.eta_det = eta_chip, eta_det
        self.Lambda_0 = Lambda_0
        self.omega_q = omega_q_hz
        self.f_cav = f_cav_hz

        self.F_max = purcell_max(Q_loaded, V_rel)
        self.F_C = 1.0 + (self.F_max - 1.0) * XI_POL * xi_pos
        # spin-flipping partner (C2) is detuned by the qubit frequency
        LC2 = lorentzian_suppression(omega_q_hz, f_cav_hz, Q_loaded)
        self.F_C2 = 1.0 + (self.F_max - 1.0) * XI_POL * xi_pos * LC2

        # lifetime-reduction factor (Lee et al. Eq. 3 with QE)
        self.zeta = 1.0 + ETA_Q * ETA_DW * ETA_BR * (self.F_C - 1.0)
        self.gamma_cav = GAMMA_0 * self.zeta
        self.tau_cav = 1.0 / self.gamma_cav

        # cyclicity: conserving rate x zeta ; flipping rate x zeta_flip
        zeta_flip = 1.0 + ETA_Q * ETA_DW * ETA_BR * (self.F_C2 - 1.0)
        self.Lambda_cav = Lambda_0 * self.zeta / zeta_flip

        # fraction of decays that put a photon in the cavity mode
        self.beta = ETA_Q * ETA_DW * ETA_BR * self.F_C / self.zeta
        # total detection efficiency per emission event
        self.eta = (self.beta * eta_wg * T_taper * eta_chip * eta_det)

    def summary(self):
        return {k: float(v) for k, v in dict(
            Q=self.Q, V=self.V, F_max=self.F_max, F_C=self.F_C,
            F_C2=self.F_C2, zeta=self.zeta, tau_cav_ns=self.tau_cav * 1e9,
            Lambda_cav=self.Lambda_cav, beta=self.beta,
            eta_wg=self.eta_wg, T_taper=self.T_taper,
            eta_chip=self.eta_chip, eta_det=self.eta_det,
            eta=self.eta).items()}


# ------------------------------------------------- readout photon model

def polarization_rate(gamma, Lambda, s=10.0, delta_over_gamma=0.0):
    """Gamma_p = R/(1+Lambda) with R = (gamma/2) s/(1+s+(2 delta/gamma)^2)
    (Rosenthal Eq. A8), s = p/psat."""
    R = 0.5 * gamma * s / (1.0 + s + (2 * delta_over_gamma) ** 2)
    return R / (1.0 + Lambda)


def readout_counts(eta, Lambda, gamma, tau, s=10.0, noise_rate=0.0):
    """(nb, nd) for a readout window tau (Rosenthal Eqs. A13-A15)."""
    Gp = polarization_rate(gamma, Lambda, s)
    nd = noise_rate * tau
    nb = nd + eta * (Lambda + 1.0) * (1.0 - np.exp(-Gp * tau))
    return nb, nd


def fidelity(nb, nd, f0=1.0):
    """Single-shot fidelity, threshold N_r = 1 (Rosenthal Eq. A20)."""
    return 0.5 + 0.5 * f0 * (np.exp(-nd) - np.exp(-nb))


def geometric_pmf_cdf(mean, kmax):
    """Detected-photon distribution for flip-terminated emission.

    Competition per emission event: a detected spin-conserving photon
    (prob p = eta*Lambda/(eta*Lambda + 1)) versus a spin flip.  The
    detected count K is geometric, P(K = k) = p^k (1 - p), with mean
    eta*Lambda.  Valid for readout windows tau >~ 3/Gamma_p (emission
    terminated by the flip, not by the window)."""
    p = mean / (mean + 1.0)
    ks = np.arange(kmax + 1)
    pmf = (1 - p) * p ** ks
    return pmf, np.cumsum(pmf)


def poisson_cdf(mean, kmax):
    ks = np.arange(kmax + 1)
    logp = ks * np.log(np.maximum(mean, 1e-300)) - mean - \
        np.array([np.sum(np.log(np.arange(1, k + 1))) for k in ks])
    pmf = np.exp(logp)
    return pmf, np.cumsum(pmf)


def _bright_ge(eta, R, Gp, tau, kmax, nt=41):
    """P(K >= Nr) for Nr = 1..kmax-1 for the exact bright record: photon
    detections are a Poisson process of rate eta*R while the spin is
    bright; the spin flips (record ends) at rate Gp.

    P(K = k) = int_0^tau Gp e^{-Gp t} Po(k; eta R t) dt
               + e^{-Gp tau} Po(k; eta R tau).
    Returns array P_ge[Nr] (index 0 unused)."""
    ts = np.linspace(1e-12, tau, nt)
    wts = np.gradient(ts) * Gp * np.exp(-Gp * ts)
    ks = np.arange(kmax)
    lgf = np.concatenate(([0.0], np.cumsum(np.log(np.arange(1, kmax)))))
    mus = eta * R * ts
    logp = (ks[:, None] * np.log(mus[None, :]) - mus[None, :]
            - lgf[:, None])
    pmf_t = np.exp(logp)                     # (kmax, nt)
    pmf = pmf_t @ wts + np.exp(-Gp * tau) * pmf_t[:, -1]
    pmf /= pmf.sum()
    cdf = np.cumsum(pmf)
    P_ge = np.empty(kmax)
    P_ge[0] = 1.0
    P_ge[1:] = 1.0 - cdf[:-1]
    return P_ge


def fidelity_exact(eta, Lambda0, gamma_cav, leak, tau,
                   noise_rate=1e3, s=1.0, f0=1.0, nt=9):
    """Single-shot fidelity, optimal integer threshold, exact bright
    statistics (Poisson process terminated by the spin flip) and dark
    statistics including background, direct off-resonant scattering and
    dark-state switch-on (off-resonant excitation followed by a spin
    flip into the bright manifold at rate R_sw = leak*R*Lambda/(1+Lambda)).
    Valid in both the window-limited (Gp tau << 1) and flip-terminated
    (Gp tau >> 1) regimes.  Returns (Fr, Nr, m_b, n_d)."""
    R = 0.5 * gamma_cav * s / (1.0 + s)
    Gp = R / (1.0 + Lambda0)
    m_b = eta * (Lambda0 + 1.0) * (1.0 - np.exp(-Gp * tau))
    m_win = eta * R * tau
    kmax = int(min(min(m_b, m_win) + 8 * np.sqrt(min(m_b, m_win)) + 50,
                   6000))
    P_ge_b = _bright_ge(eta, R, Gp, tau, kmax)
    n_d = eta * leak * R * tau + noise_rate * tau
    _, cdf_d = poisson_cdf(n_d, kmax)
    # dark switch-on: switch at t, then a bright record of length tau-t
    R_sw = leak * R * Lambda0 / (Lambda0 + 1.0)
    ts = np.linspace(1e-12, tau, nt)
    wts = np.gradient(ts) * R_sw * np.exp(-R_sw * ts)
    P_sw = np.zeros(kmax)
    for t, w in zip(ts, wts):
        P_sw += w * _bright_ge(eta, R, Gp, tau - t + 1e-15, kmax, nt=15)
    Nrs = np.arange(1, kmax)
    err_dark = np.minimum(1.0, (1.0 - cdf_d[Nrs - 1]) + P_sw[Nrs])
    err_bright = 1.0 - P_ge_b[Nrs]
    Frs = 1.0 - 0.5 * (err_dark + err_bright)
    i = int(np.argmax(Frs))
    return f0 * float(Frs[i]) + (1 - f0) * 0.5, int(Nrs[i]), m_b, n_d


def fidelity_threshold(eta, Lambda0, gamma_cav, leak, tau,
                       noise_rate=1e3, s=1.0, f0=1.0):
    """Single-shot fidelity with optimal integer threshold N_r.

    Bright (|dn>): the record is flip-terminated -- detected count
      K ~ Geometric with mean m_b = eta (Lambda0+1)(1 - exp(-Gp tau)),
      P(K >= Nr) = p^Nr with p = m_b/(m_b+1).
    Dark (|up>): three contributions --
      (a) background Poisson(noise_rate * tau);
      (b) direct off-resonant scattering photons, Poisson mean
          eta * leak * R * tau;
      (c) SWITCH-ON: an off-resonant excitation of the dark spin decays
          spin-flippingly into the bright manifold with probability
          Lambda/(Lambda+1), after which the emitter turns bright for
          the remainder of the window.  Switch rate R_sw = leak * R.
          A switch at time t contributes a bright-like record with mean
          m(t) = eta (Lambda0+1)(1 - exp(-Gp (tau-t))).
    Here R = (gamma_cav/2) s/(1+s) and leak = R_dark/R (power-broadening
    included by the caller).  Returns (Fr, Nr_opt, m_b, n_d_direct)."""
    R = 0.5 * gamma_cav * s / (1.0 + s)
    Gp = R / (1.0 + Lambda0)
    m_b = eta * (Lambda0 + 1.0) * (1.0 - np.exp(-Gp * tau))
    p = m_b / (m_b + 1.0)
    R_sw = leak * R * Lambda0 / (Lambda0 + 1.0)
    n_d = eta * leak * R * tau + noise_rate * tau       # direct + bg
    kmax = int(max(30, 5 * n_d + 60))
    _, cdf_d = poisson_cdf(n_d, kmax)
    # switch-on: quadrature over switch time (exponential arrival)
    ts = np.linspace(0.0, tau, 41)
    w = np.gradient(ts) * R_sw * np.exp(-R_sw * ts)     # arrival density
    m_t = eta * (Lambda0 + 1.0) * (1.0 - np.exp(-Gp * (tau - ts)))
    p_t = m_t / (m_t + 1.0)
    Nrs = np.arange(1, kmax)
    err_dark_po = 1.0 - cdf_d[Nrs - 1]
    err_dark_sw = (w[None, :] * p_t[None, :] ** Nrs[:, None]).sum(axis=1)
    err_dark = np.minimum(1.0, err_dark_po + err_dark_sw)
    err_bright = 1.0 - p ** Nrs
    Frs = 1.0 - 0.5 * (err_dark + err_bright)
    i = int(np.argmax(Frs))
    return f0 * float(Frs[i]) + (1 - f0) * 0.5, int(Nrs[i]), m_b, n_d


def validate_confocal():
    """Reproduce the measured confocal operating point of 2403.13110:
    eta ~ 0.2%, Lambda = 2244, tau = 50 us -> nb ~ 4, Fr = 87.4%
    (f0 inferred from state-preparation errors)."""
    eta, Lam, tau = 0.002, 2244.0, 50e-6
    nb, nd = readout_counts(eta, Lam, GAMMA_0, tau, s=10.0,
                            noise_rate=0.2 / 50e-6)
    # nd fixed to the reported <= 0.2 counts via the noise rate above
    f0 = 0.935
    Fr = fidelity(nb, nd, f0)
    return dict(nb=float(nb), nd=float(nd), Fr=float(Fr),
                target_nb=4.0, target_Fr=0.874)


if __name__ == "__main__":
    print("confocal validation:", json.dumps(validate_confocal(), indent=1))

    # design point (numbers refined by cavity/EME stages)
    ci = CavityInterface(Q_loaded=1.0e5, V_rel=0.73, eta_wg=0.85)
    print(json.dumps(ci.summary(), indent=1))
    for tau_us in (0.1, 0.5, 1, 5):
        nb, nd = readout_counts(ci.eta, ci.Lambda_cav, ci.gamma_cav,
                                tau_us * 1e-6, noise_rate=1e4)
        print(f"tau={tau_us:5.1f} us: nb={nb:8.1f} nd={nd:.3f} "
              f"Fr={fidelity(nb, nd, 0.99):.6f}")
