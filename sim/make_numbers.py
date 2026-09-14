"""Collect every number quoted in the manuscript from the results files
and write LaTeX macros.  The macro block is written to paper/macros.tex
and also spliced into paper/manuscript.tex and supplement/supplement.tex
between the AUTO-NUMBERS markers if those files exist."""

import json
import pathlib
import re

import numpy as np

SIM = pathlib.Path(__file__).resolve().parent
ROOT = SIM.parent
RES = SIM / "results"

import sys
sys.path.insert(0, str(SIM))
import materials
import readout as ro
import interface as itf

r = {}

# materials
r["nDia"] = f"{float(materials.diamond.n(0.619)):.3f}"
r["nLNo"] = f"{float(materials.linbo3_o.n(0.619)):.2f}"
r["nLNe"] = f"{float(materials.linbo3_e.n(0.619)):.2f}"
r["nSiO"] = f"{float(materials.silica.n(0.619)):.3f}"

# waveguide / taper
eme = json.load(open(RES / "taper_eme.json"))
i6 = eme["lengths"].index(6.0)
r["TaperT"] = f"{100*eme['T'][i6]:.1f}"
r["TaperTthree"] = f"{100*eme['T'][eme['lengths'].index(3.0)]:.1f}"
r["TaperLen"] = "6"
r["WgWidth"] = "280"
r["LNWidth"] = "600"

# cavity (production, gmax 2.5)
prod = json.load(open(RES / "cavity_production.json"))
r["CavA"] = f"{prod['a_nm']:.1f}"
r["CavLam"] = f"{prod['lam_nm']:.1f}"
r["CavQi"] = f"{prod['Q']/1e6:.1f}"
r["MirrorD"] = f"{prod['mirror_decay_D']:.2f}"
r["MirrorDdB"] = f"{10*np.log10(1/prod['mirror_decay_D']):.1f}"
r["CavV"] = "0.68"
r["CavVspread"] = "0.65--0.78"
r["HoleR"] = f"{0.30*prod['a_nm']:.0f}"

# convergence spread
import glob
qvals = [json.load(open(f))["Q"] for f in glob.glob(str(RES/"cavity_mode_a193*.json"))
         if json.load(open(f))["N_mirror"] == 10]
r["CavQspread"] = f"{100*(max(qvals)/min(qvals)-1):.0f}"

# tolerances
tol = json.load(open(RES / "tolerance.json"))
lam0 = [t for t in tol if t["tag"] == "nominal_g2"][0]["lam_nm"]
sens = {t["tag"]: t["lam_nm"] - lam0 for t in tol}
r["TolR"] = f"{abs(sens['r_plus'] - sens['r_minus'])/2/2:.1f}"   # nm per nm radius
r["TolW"] = f"{abs(sens['w_plus'] - sens['w_minus'])/2/10:.2f}"
r["QSuspended"] = f"{[t for t in tol if t['tag']=='suspended'][0]['Q']/1e6:.1f}"

# spin model
snvr = json.load(open(RES / "snv_results.json"))
r["DeltaG"] = f"{snvr['delta_g_GHz']:.0f}"
r["FqModel"] = f"{snvr['fq_147_125mT']:.2f}"
r["CycMinModel"] = f"{snvr['cyc_53_180mT']:.1f}"
r["DevStrain"] = f"{snvr['device_strain_ratio']:.2f}"

# design point
dp = itf.design(500.0)
r["QL"] = "500"
r["KappaTHz"] = f"{dp['kappa_GHz']/1e3:.1f}"
r["FC"] = f"{dp['F_C']:.0f}"
r["FD"] = f"{dp['F_D']:.1f}"
r["Zeta"] = f"{dp['zeta']:.1f}"
r["TauCav"] = f"{dp['tau_cav_ns']:.2f}"
r["GammaCav"] = f"{dp['gamma_cav_GHz']*1e3:.0f}"     # MHz
r["Beta"] = f"{100*dp['beta']:.0f}"
r["EtaWg"] = f"{100*dp['eta_wg']:.2f}"
r["EtaTot"] = f"{100*dp['eta']:.0f}"
r["gGHz"] = f"{dp['g_GHz']:.1f}"

# readout at design point
itfj = json.load(open(RES / "interface.json"))
ro_ = itfj["readout"]
r["FrMax"] = f"{100*ro_['max_cyclicity']['Fr']:.1f}"
r["FrMaxTau"] = f"{ro_['max_cyclicity']['tau_us']:.2f}"
r["FrMaxNb"] = f"{ro_['max_cyclicity']['nb']:.0f}"
r["FrStrain"] = f"{100*ro_['high_strain']['Fr']:.1f}"
r["FrStrainTau"] = f"{1e3*ro_['high_strain']['tau_us']:.0f}"     # ns
r["FrWorst"] = f"{100*ro_['worst_angle']['Fr']:.1f}"
r["QLoptHundred"] = f"{itfj['QL_opt']['100.0']['QL']:.0f}"
r["FrOptHundred"] = f"{100*itfj['QL_opt']['100.0']['Fr']:.1f}"
r["FrOptMax"] = f"{100*itfj['QL_opt']['2244.0']['Fr']:.1f}"

# comparison
frl = json.load(open(RES / "fr_vs_lambda.json"))
r["EtaConfocal"] = "0.2--0.4"

# verification
ver = json.load(open(RES / "verification.json"))
r["NChecks"] = str(ver["_summary"]["total"] - 0)
r["NChecksPassed"] = str(ver["_summary"]["passed"])

lines = ["% ---- auto-generated result macros (source: sim/make_numbers.py) ----"]
for k, v in r.items():
    lines.append(f"\\newcommand{{\\n{k}}}{{{v}}}")
lines.append("% ---- end of result macros ----")
block = "\n".join(lines)

(ROOT / "paper").mkdir(exist_ok=True)
(ROOT / "paper" / "macros.tex").write_text(block + "\n")
print(block)

for tex in (ROOT / "paper" / "manuscript.tex",
            ROOT / "supplement" / "supplement.tex"):
    if tex.exists():
        s = tex.read_text()
        s2 = re.sub(r"% ---- auto-generated result macros.*?% ---- end of result macros ----",
                    lambda m: block, s, flags=re.S)
        tex.write_text(s2)
        print("spliced into", tex)
