"""Production SnV Hamiltonian results: validation + strain trade-off data."""
import json, pathlib
import numpy as np
import snv

RES = pathlib.Path("results"); RES.mkdir(exist_ok=True)
prm = snv.PARAMS
th, ph = prm['theta'], prm['phi']
mu = np.array([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)])
e1 = np.cross(mu, [0,0,1.0]); e1 /= np.linalg.norm(e1)
e2 = np.cross(mu, e1)

def b_misaligned(alpha_deg, B):
    a = np.deg2rad(alpha_deg)
    return B*(np.cos(a)*mu + np.sin(a)*e1)

out = {}

# 1) validation numbers
out['delta_g_GHz'] = float(np.sqrt(prm['lam_g']**2 + 4*prm['ups_g']**2))
out['delta_e_GHz'] = float(np.sqrt(prm['lam_e']**2 + 4*prm['ups_e']**2))
out['fq_147_125mT'] = snv.qubit_frequency(snv.b_lab(147.0, 0.125))
out['cyc_53_180mT'] = snv.cyclicity(snv.b_lab(53.0, 0.180))
out['cyc_147_180mT'] = snv.cyclicity(snv.b_lab(147.0, 0.180))

# 2) cyclicity + qubit freq vs zeta at 180 mT (validation curve)
zs = np.arange(0.0, 360.0, 2.0)
out['zeta_deg'] = list(zs)
out['cyc_vs_zeta'] = [snv.cyclicity(snv.b_lab(z, 0.180)) for z in zs]
out['fq_vs_zeta_125mT'] = [snv.qubit_frequency(snv.b_lab(z, 0.125)) for z in zs]

# 3) cyclicity vs misalignment angle alpha (at 180 mT)
alphas = np.concatenate([np.arange(0.5, 5, 0.5), np.arange(5, 46, 1.0)])
out['alpha_deg'] = list(alphas)
out['cyc_vs_alpha'] = [snv.cyclicity(b_misaligned(a, 0.180)) for a in alphas]

# 4) strain trade-off: Lambda and Rabi vs 2*Ups_g/lambda_g at alpha=9 deg
#    (drive 0.6 mT perpendicular to the dipole; B = 180 mT)
strain_ratio = np.logspace(np.log2(0.02), np.log2(8.0), 40, base=2)
scales = strain_ratio * prm['lam_g'] / (2*prm['ups_g'])
Bv = b_misaligned(9.0, 0.180)
bdrive = 6e-4 * e1
out['strain_2UgOverLam'] = list(strain_ratio)
out['cyc_vs_strain'] = [snv.cyclicity(Bv, ups_scale=s) for s in scales]
out['rabi_MHz_vs_strain'] = [snv.rabi_rate(Bv, bdrive, ups_scale=s)*1e3 for s in scales]
out['fq_GHz_vs_strain'] = [snv.qubit_frequency(Bv, ups_scale=s) for s in scales]
# device point
dev_ratio = 2*prm['ups_g']/prm['lam_g']
out['device_strain_ratio'] = float(dev_ratio)
out['device_cyc'] = snv.cyclicity(Bv)
out['device_rabi_MHz'] = snv.rabi_rate(Bv, bdrive)*1e3

json.dump({k: (v if not isinstance(v, list) else [float(x) for x in v])
           for k, v in out.items()}, open(RES/'snv_results.json','w'), indent=1)
print("delta_g", round(out['delta_g_GHz'],2), "fq", round(out['fq_147_125mT'],3),
      "cyc53", round(out['cyc_53_180mT'],1), "cyc@9deg",
      round(snv.cyclicity(b_misaligned(9.0,0.180)),0))
print("device: ratio", round(dev_ratio,3), "cyc", round(out['device_cyc'],0),
      "rabi", round(out['device_rabi_MHz'],2), "MHz")
print("strain sweep endpoints: cyc", round(out['cyc_vs_strain'][0],0), "->",
      round(out['cyc_vs_strain'][-1],1), "; rabi", round(out['rabi_MHz_vs_strain'][0],3),
      "->", round(out['rabi_MHz_vs_strain'][-1],2), "MHz")
