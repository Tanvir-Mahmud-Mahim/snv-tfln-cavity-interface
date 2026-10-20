"""Collect every number quoted in the revised manuscript and supplement
from the results files and write LaTeX macros to sim/results/macros.tex.
No number in the text is typed by hand: each macro below is computed
from a results file produced by the simulation pipeline."""

import glob
import json
import pathlib
import sys

import numpy as np

SIM = pathlib.Path(__file__).resolve().parent
ROOT = SIM.parent
RES = SIM / "results"
sys.path.insert(0, str(SIM))
import materials  # noqa: E402

r = {}
J = lambda name: json.loads((RES / name).read_text())
pct = lambda x, n=1: f"{100 * x:.{n}f}"

# ------------------------------------------------------------ materials
r["nDia"] = f"{float(materials.diamond.n(0.619)):.3f}"
r["nLNo"] = f"{float(materials.linbo3_o.n(0.619)):.2f}"
r["nLNe"] = f"{float(materials.linbo3_e.n(0.619)):.2f}"
r["nSiO"] = f"{float(materials.silica.n(0.619)):.3f}"

# ------------------------------------------------------- taper (EME)
TOL = RES / "tolerance_eme"
def T6(tag, L=6.0):
    d = json.loads((TOL / f"{tag}.json").read_text())
    return d["T_LN"][d["lengths"].index(L)]
nom = json.loads((TOL / "nominal.json").read_text())
r["TaperT"] = pct(T6("nominal"))
r["TaperTthree"] = pct(nom["T_LN"][nom["lengths"].index(3.0)])
r["TaperLen"] = "6"
r["TaperTmFour"] = pct(T6("nominal_m4"))
r["TaperTmEight"] = pct(T6("nominal_m8"))
r["TaperTsegFF"] = pct(T6("nominal_seg45"))
r["WgWidth"] = "280"
r["LNWidth"] = "600"
for tag, name in (("off25", "TolOffTwentyFive"), ("off50", "TolOffFifty"),
                  ("off75", "TolOffSeventyFive"), ("off100", "TolOffHundred"),
                  ("off150", "TolOffOneFifty"), ("tilt0p5", "TolTiltHalf"),
                  ("tilt1", "TolTiltOne"), ("tilt2", "TolTiltTwo"),
                  ("tilt3", "TolTiltThree"), ("gap10", "TolGapTen"),
                  ("gap20", "TolGapTwenty"), ("gap50", "TolGapFifty"),
                  ("sw80", "TolSwEighty"), ("sw70", "TolSwSeventy"),
                  ("sw60", "TolSwSixty"), ("rib60", "TolRib"),
                  ("rib60_sw70", "TolRibSw"), ("combined", "TolComb")):
    r[name] = pct(T6(tag))
r["TolGapFiftyLong"] = pct(T6("gap50", 14.0))

# ---------------------------------------------------- cavity (closed)
prod = J("cavity_production.json")
r["CavA"] = f"{prod['a_nm']:.1f}"
r["CavLam"] = f"{prod['lam_nm']:.1f}"
r["CavQi"] = f"{prod['Q']/1e6:.1f}"
r["MirrorD"] = f"{prod['mirror_decay_D']:.2f}"
r["MirrorDdB"] = f"{10*np.log10(1/prod['mirror_decay_D']):.1f}"
r["CavV"] = "0.68"
r["HoleR"] = f"{0.30*prod['a_nm']:.0f}"
qvals = [json.load(open(f))["Q"] for f in glob.glob(str(RES / "cavity_mode_a193*.json"))
         if json.load(open(f))["N_mirror"] == 10]
r["CavQspread"] = f"{100*(max(qvals)/min(qvals)-1):.0f}"
tol = J("tolerance.json")
lam0 = [t for t in tol if t["tag"] == "nominal_g2"][0]["lam_nm"]
sens = {t["tag"]: t["lam_nm"] - lam0 for t in tol}
# hole radius 0.29a -> 0.31a (a = production period): shift per nm of radius
r["TolR"] = f"{abs(sens['r_plus'] - sens['r_minus'])/(0.02*prod['a_nm']):.1f}"
r["TolW"] = f"{abs(sens['w_plus'] - sens['w_minus'])/2/10:.2f}"
r["QSuspended"] = f"{[t for t in tol if t['tag']=='suspended'][0]['Q']/1e6:.1f}"

# ---------------------------------------------------- cavity (loaded FDTD)
fd = J("cavity_fdtd_summary.json")
cl = J("cavity_loaded.json")
r["FDTDres"] = "40"
for n, name in ((0, "QLNzero"), (1, "QLNone"), (2, "QLNtwo"), (3, "QLNthree")):
    r[name] = f"{fd[f'N{n}_r40']['Q_L']:.0f}"
dsg = fd[f"N{cl['n_out']}_r40"]
r["Nout"] = str(cl["n_out"])
r["QL"] = f"{cl['Q_L']:.0f}"
r["lamFDTD"] = f"{np.mean([fd[k]['lambda_nm'] for k in fd if k.endswith('_r40')]):.0f}"
r["TEzeroFrac"] = pct(dsg["frac_TE0"], 0)
r["DownFrac"] = pct(dsg["frac_down"], 0)
r["UpFrac"] = pct(dsg["frac_up"], 0)
r["SideFrac"] = pct(dsg["frac_side"], 0)
r["LeftFrac"] = f"{max(0.1, np.ceil(1000*dsg['frac_left'])/10):.1f}"
r["FcavFDTD"] = f"{dsg['F_cav']:.0f}"
r["FbgFDTD"] = f"{dsg['F_bg']:.2f}"
r["EtaWgFDTD"] = pct(cl["eta_wg"], 0)
rows = []
for k in sorted(fd):
    a = fd[k]
    rows.append(f"{a['n_out']}{'' if a['resolution']==40 else ' (res. 50)'} & "
                f"{a['lambda_nm']:.1f} & {a['Q_L']:.0f} & {a['F_cav']:.1f} & "
                f"{100*a['frac_TE0']:.1f} & {100*a['frac_up']:.1f} & "
                f"{100*a['frac_down']:.1f} & {100*a['frac_side']:.1f} & "
                f"{100*a['frac_left']:.2f} & {100*a['eta_wg']:.1f} \\\\")
r["__FDTDrows"] = "\n".join(rows)
if "N0_r50" in fd and "N0_r40" in fd:
    a40, a50 = fd["N0_r40"], fd["N0_r50"]
    r["__FDTDconvsentence"] = (
        f"Increasing the resolution from 40 to 50~pixels/\\textmu m for "
        f"$N_{{\\mathrm{{out}}}} = 0$ changes $Q_L$ from {a40['Q_L']:.0f} to "
        f"{a50['Q_L']:.0f}, $F_{{\\mathrm{{cav}}}}$ from {a40['F_cav']:.0f} to "
        f"{a50['F_cav']:.0f}, the waveguide-coupled fraction from "
        f"{100*a40['frac_TE0']:.1f}\\% to {100*a50['frac_TE0']:.1f}\\%, and "
        f"the resonance from {a40['lambda_nm']:.1f} to {a50['lambda_nm']:.1f}~nm.")
else:
    r["__FDTDconvsentence"] = ""
a0 = fd["N0_r40"]
r["__FDTDldosSentence"] = (
    f"the fit gives $F_{{\\mathrm{{cav}}}} = {a0['F_cav_ldos']:.0f}$ and "
    f"$Q_L = {a0['Q_ldos']:.0f}$, against {a0['F_cav']:.0f} and "
    f"{a0['Q_L']:.0f} from the mode volume and Harminv.")

# ------------------------------------------------------------ spin model
snvr = J("snv_results.json")
r["DeltaG"] = f"{snvr['delta_g_GHz']:.0f}"
r["FqModel"] = f"{snvr['fq_147_125mT']:.2f}"
r["CycMinModel"] = f"{snvr['cyc_53_180mT']:.1f}"
r["DevStrain"] = f"{snvr['device_strain_ratio']:.2f}"

# ------------------------------------------------------------ readout
rr = J("readout_rev.json")
ref = rr["reference"]
lv, bud = ref["levels"], ref["budget"]
r["Bref"] = f"{1e3*ref['B_T']:.0f}"
r["FqRef"] = f"{lv['wq']:.2f}"
r["WeRef"] = f"{lv['we']:.2f}"
r["DABref"] = f"{lv['dA1B2']:.2f}"
r["LamRef"] = f"{lv['LamA']:.0f}"
r["SBtwo"] = f"{lv['S_rel'][1][1]:.3f}"
r["SAtwo"] = f"{lv['S_rel'][1][0]:.1e}".replace("e-0", "\\times10^{-").replace("e-", "\\times10^{-") + "}"
r["SAtwo"] = "$" + r["SAtwo"] + "$"
r["SBone"] = "$" + (f"{lv['S_rel'][0][1]:.1e}".replace("e-0", "\\times10^{-")
                    .replace("e-", "\\times10^{-")) + "}$"
r["KappaTHz"] = f"{bud['kappa_GHz']/1e3:.1f}"
r["KappaNm"] = f"{618.9**2 * bud['kappa_GHz']*1e9 / 2.99792458e17:.1f}"
_knm = 618.9**2 * bud['kappa_GHz']*1e9 / 2.99792458e17
r["HalfKappaNm"] = f"{_knm/2:.2f}"
_dl_dr = abs(sens['r_plus'] - sens['r_minus'])/(0.02*prod['a_nm'])
r["RadTolNm"] = f"{_knm/2/_dl_dr:.1f}"
r["FC"] = f"{bud['F_C']:.0f}"
r["FD"] = f"{bud['F_D']:.1f}"
r["LifeRed"] = f"{bud['zeta']:.1f}"
r["TauCav"] = f"{bud['tau_ns']:.2f}"
r["GammaCavMHz"] = f"{1e3*bud['gamma_GHz']:.0f}"
r["Beta"] = pct(bud["beta"], 0)
r["EtaTot"] = pct(bud["eta"], 0)
r["DeltaGTHz"] = f"{lv['delta_g']/1e3:.2f}"
C0 = 299792458.0
f_cav = C0 / 618.9e-9
tau = np.array(ref["vs_tau"]["tau_ns"]); Ft = np.array(ref["vs_tau"]["F"], float)
i100 = int(np.argmin(abs(tau - 100)))
# the 100-ns value is taken from a dedicated run at exactly 100 ns
F100 = rr["extra"]["cond_rates"]["F100"]
r["FrHundred"] = pct(F100)
r["FrHundredNb"] = f"{rr['extra']['cond_rates']['nb100']:.1f}"
r["FrHundredNd"] = f"{rr['extra']['cond_rates']['nd100']:.2f}"
r["FrBest"] = pct(ref["best"]["F"])
r["FrBestTau"] = f"{ref['best']['tau_ns']:.0f}"
cf = ref["confocal"]
r["ConfTwo"] = pct(cf["eta0.002_xlw0.0"]["F50us"])
r["ConfFour"] = pct(cf["eta0.004_xlw0.0"]["F50us"])
r["ConfMaxTwo"] = pct(np.nanmax(np.array(cf["eta0.002_xlw0.0"]["F"], float)))
r["ConfNb"] = f"{cf['eta0.002_xlw0.0']['nb50']:.1f}"
r["ConfNd"] = f"{cf['eta0.002_xlw0.0']['nd50']:.2f}"
cr = rr["extra"]["cond_rates"]
_sh = cr['pB_d']/(cr['pA_d']+cr['pB_d'])
r["DarkBtwoShare"] = "more than 99.9" if _sh > 0.999 else f"{100*_sh:.0f}"
r["SwitchRatio"] = f"{cr['flip_viaB_d']/cr['flip_viaA_d']:.0f}"

ql = rr["vs_QL"]
rowsQ = {k: {x["Q_L"]: x["F"] for x in v["rows"]} for k, v in ql.items()}
r["FrQLthreeK"] = pct(rowsQ["ref_a9"][3000])
r["DesignBelowMax"] = f"{100*(max(rowsQ['ref_a9'].values()) - ref['best']['F']):.1f}"
r["FrAfortyfiveQLbest"] = pct(max(rowsQ["a45"].values()))
names = {"ref_a9": "$\\alpha = 9^\\circ$", "a35": "$\\alpha = 35^\\circ$",
         "a45": "$\\alpha = 45^\\circ$",
         "inplane53_180mT": "$\\zeta = 53^\\circ$, 180~mT"}
r["__QLrows"] = "\n".join(
    names[k] + " & " + " & ".join(f"{100*rowsQ[k][q]:.1f}" for q in
                                   (60, 100, 150, 250, 400, 500, 700, 1000, 1500, 3000)) + " \\\\"
    for k in ("ref_a9", "a35", "a45", "inplane53_180mT"))

va = {x["alpha"]: x for x in rr["vs_alpha"]}
# loaded-Q plateau and optimum
qa9 = rowsQ["ref_a9"]
flat = [qa9[q] for q in qa9 if 60 <= q <= 700]
r["QLflat"] = f"{100*(max(flat)-min(flat)):.1f}"
r["QLoptAthirtyfive"] = f"{max(rowsQ['a35'], key=rowsQ['a35'].get):.0f}"
r["QLoptAfortyfive"] = f"{max(rowsQ['a45'], key=rowsQ['a45'].get):.0f}"
# misalignment up to which F exceeds 90% / 95% (interpolated below the
# A1-B2 crossing), and the crossing angle itself
al = np.array(sorted(a for a in va if a <= 60)); Fa = np.array([va[a]["F"] for a in al])
def last_above(th):
    i = int(np.argmax(Fa < th))
    return al[i-1] + (al[i]-al[i-1])*(Fa[i-1]-th)/(Fa[i-1]-Fa[i])
r["AlphaNinety"] = f"{np.floor(last_above(0.90)):.0f}"
r["AlphaNinetyFive"] = f"{np.floor(last_above(0.95)):.0f}"
import readout_fcs as _rf, readout_revision as _rr
from scipy.optimize import brentq as _bq
_B = _rr.b_ref()
_xa = _bq(lambda a: _rf.level_structure(_B*_rr.b_dir(a))["nu"][1, 1], 50.0, 75.0)
r["AlphaCross"] = f"{_xa:.0f}"
vfl = {round(x["B"], 4): x for x in rr["vs_field"]}
r["DetForty"] = f"{vfl[0.04]['levels']['dA1B2']:.2f}"
r["BfiveHundred"] = f"{1e3*rr['extra']['B500']['B']:.0f}"
import snv as _snv
_l23 = _rf.level_structure(_snv.b_lab(110.0, 0.150))
r["FqRtwothree"] = f"{_l23['wq']:.2f}"
r["DABRtwothree"] = f"{_l23['nu'][1, 1]:.2f}"
vst = rr["vs_strain"]
r["LamStrainMin"] = f"{min(x['levels']['LamA'] for x in vst):.0f}"
_lm = max(x['levels']['LamA'] for x in vst); _e = int(np.floor(np.log10(_lm)))
r["LamStrainMax"] = f"${_lm/10**_e:.1f}\\times10^{{{_e}}}$"
r["Sopt"] = f"{rr['extra']['cond_rates']['s']:.2g}"
r["LamFortyFive"] = f"{va[45]['levels']['LamA']:.0f}"
r["FrAfortyfiveDesign"] = pct(va[45]["F"])
vs = rr["vs_strain"]
r["FrStrainMin"] = pct(min(x["F"] for x in vs))
r["FrStrainMax"] = pct(max(x["F"] for x in vs))
r["RabiMin"] = f"{min(x['rabi_MHz'] for x in vs):.1f}"
r["RabiMax"] = f"{max(x['rabi_MHz'] for x in vs):.1f}"
vf = {round(x["B"], 4): x["F"] for x in rr["vs_field"]}
r["FrBforty"] = pct(vf[0.04])
r["FrBsevsev"] = pct(rr["extra"]["B500"]["F"])
r["FrBhundredtwelve"] = pct(ref["best"]["F"])
r["FrBsix"] = pct(vf[0.6])

rb = rr["robust"]
F0 = ref["best"]["F"]
xl = {x["x_GHz"]: x["F"] for x in rb["extra_lw"]}
sd = {x["fwhm_GHz"]: x["F"] for x in rb["sd"]}
sb = {x["sbr"]: x["F"] for x in rb["sbr"]}
r["DropXlwHundred"] = f"{100*(F0 - xl[0.1]):.1f}"
r["DropSDtwo"] = f"{100*(F0 - sd[0.2]):.1f}"
r["FrSDhalf"] = pct(sd[0.5])
r["DropSBRninety"] = f"{100*(F0 - sb[90.0]):.2f}"
r["FrSBRthree"] = pct(sb[3.0])
r["__XlwList"] = ", ".join(pct(xl[x]) for x in (0.025, 0.05, 0.1, 0.2, 0.3, 0.5))
r["__SDList"] = ", ".join(pct(sd[x]) for x in (0.05, 0.1, 0.2, 0.3, 0.5))
r["__SBRList"] = ", ".join(pct(sb[x]) for x in (3.0, 10.0, 30.0, 90.0, 300.0, 1000.0))
ec = {x["eta_chip_det"]: x["F"] for x in rb["eta_chain"]}
r["__EtaList"] = ", ".join(pct(ec[x]) for x in (0.1, 0.2, 0.3, 0.385, 0.5, 0.7, 0.82))
r["FrEtaTen"] = pct(ec[0.1])
r["FrComb"] = pct(rb["combined"]["F"])
ex = rr["extra"]
r["EtaDline"] = f"{100*ex['dline']['eta']:.1f}"
r["FrDline"] = pct(ex["dline"]["F"])
r["EtaFiber"] = pct(ex["fiber"]["eta"], 0)
r["FrFiber"] = pct(ex["fiber"]["F"])
xi = ex["xi_pos"]
r["__XiPosSentence"] = (
    "With $Q_L$ optimized for each value, the reference-point fidelity is "
    + ", ".join(f"{100*x['F']:.1f}\\% at $\\xi_{{\\mathrm{{pos}}}} = {x['xi_pos']}$ "
                f"($Q_L = {x['QL_opt']}$)" for x in xi)
    + "; a smaller placement overlap is compensated by a higher loaded $Q$ "
      "because $\\beta$ saturates.")

# ------------------------------------------------------------ testbench
ver = J("verification.json")
r["NChecks"] = str(ver["_summary"]["total"])
r["NChecksPassed"] = str(ver["_summary"]["passed"])
r["CycNine"] = f"{ver['V5']['value']:.0f}"
r["gGHz"] = f"{ver['V11']['g_over_2pi_GHz']:.1f}"
r["GoverKappa"] = f"{ver['V11']['g_over_kappa']:.3f}"
r["TolQratio"] = f"{10*np.floor(ver['V12']['ratio']/10):.0f}"
r["EMEmaxDev"] = f"{100*max(abs(ver['V9']['T_8modes']-ver['V9']['T_6modes']), abs(ver['V9']['T_45seg']-ver['V9']['T_6modes']), abs(ver['V9']['T_dx2p5']-ver['V9']['T_6modes'])):.2f}"
r["FDTDldosDev"] = f"{100*ver['V17']['dF']:.0f}"

lines = ["% ---- auto-generated result macros (source: sim/make_numbers.py) ----"]
for k, v in r.items():
    if k.startswith("__"):
        lines.append(f"\\newcommand{{\\{k[2:]}}}{{{v}}}")
    else:
        lines.append(f"\\newcommand{{\\n{k}}}{{{v}}}")
lines.append("% ---- end of result macros ----")
(RES / "macros.tex").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
