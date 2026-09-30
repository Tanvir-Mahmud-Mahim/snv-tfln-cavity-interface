"""Readout fidelity versus bare cyclicity Lambda0 (Fig. 1(c), Fig. 5(b)).

cavity:        interface.readout_point() at the design point Q_L = 500
               (same optimisation as interface.py); cavity_tau_us is the
               chosen readout window.
confocal_002,  readout.fidelity_exact() for a confocal set-up without cavity
confocal_004:  (detection efficiency 0.2% and 0.4%), with the settings of
               check V10 in checks.py (bulk decay rate, s = 2, 35.4 MHz
               linewidth, 0.15 background counts per 50 us) and f0 = 0.99;
               the better of a 10 us and a 30 us window.
The archived file (Zenodo 10.5281/zenodo.22758086) does not record these
settings; they were recovered by comparison with it and reproduce it
exactly.  Writes results/fr_vs_lambda.json."""
import json

import numpy as np

import interface as itf
import readout as ro

LAMBDA0 = np.logspace(np.log10(3.0), np.log10(5000.0), 16)
S_CONF = 2.0
LEAK_CONF = (1 + S_CONF) / (1 + S_CONF + (2 * 3.677e9 / 35.4e6) ** 2)
TAUS_CONF = (10e-6, 30e-6)

dp = itf.design(500.0)
out = {"Lambda0": [float(x) for x in LAMBDA0], "cavity": [],
       "cavity_tau_us": [], "confocal_002": [], "confocal_004": []}
for lam0 in LAMBDA0:
    nb, nd, Fr, Gp, topt, Nr, sopt = itf.readout_point(dp, float(lam0))
    out["cavity"].append(float(Fr))
    out["cavity_tau_us"].append(float(topt * 1e6))
    for eta_c, key in ((0.002, "confocal_002"), (0.004, "confocal_004")):
        out[key].append(max(
            ro.fidelity_exact(eta_c, float(lam0), ro.GAMMA_0, LEAK_CONF, t,
                              noise_rate=0.15 / 50e-6, s=S_CONF, f0=0.99)[0]
            for t in TAUS_CONF))
    print(f"L0={lam0:7.1f}  cavity Fr={out['cavity'][-1]:.4f} "
          f"(tau={out['cavity_tau_us'][-1]:.3f} us)  confocal Fr="
          f"{out['confocal_002'][-1]:.4f}-{out['confocal_004'][-1]:.4f}",
          flush=True)
(itf.RESULTS / "fr_vs_lambda.json").write_text(json.dumps(out, indent=1))
print("saved fr_vs_lambda.json")
