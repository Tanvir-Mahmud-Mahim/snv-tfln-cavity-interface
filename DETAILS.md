# Fast Single-Shot Readout of Tin-Vacancy Spins: details

This file has the full notes for version 2.0.1 (revised manuscript). For
the overview, installation and quick start, see [README.md](README.md).

## Contents

- [1. What the code computes](#1-what-the-code-computes)
- [2. The full file tree](#2-the-full-file-tree)
- [3. Fonts](#fonts)
- [4. Files the figure and number scripts expect](#4-files-the-figure-and-number-scripts-expect)
- [5. The scripts, step by step](#5-the-scripts-step-by-step)
- [6. Which script makes which figure](#6-which-script-makes-which-figure)
- [7. The Python modules](#7-the-python-modules)
- [8. Where the numbers come from](#8-where-the-numbers-come-from)
- [9. Built-in checks](#9-built-in-checks)
- [10. Notes on the calculations](#10-notes-on-the-calculations)
- [11. Version history](#11-version-history)

---

## 1. What the code computes

The code follows the light from the emitter to the detector:

- refractive indices of diamond, lithium niobate and SiO<sub>2</sub> at
  619 nm;
- the modes of the diamond beam, the TFLN ridge and the hybrid stack
  (finite elements), and the taper transfer from diamond to TFLN
  (eigenmode expansion, EME), including lateral offset, angular
  misalignment, an air gap, sloped sidewalls and a partially etched ridge;
- the closed nanobeam cavity (guided-mode expansion, GME): resonance,
  intrinsic Q, mode volume, fabrication tolerances;
- the loaded one-sided cavity (3D FDTD with Meep): loaded Q, Purcell
  factor, and the fractions of the emission that reach the output
  waveguide, the substrate, air and the closed mirror;
- the SnV<sup>-</sup> spin Hamiltonian, checked against published
  measurements (energy levels, qubit frequency, cyclicity versus field
  angle);
- the readout: a four-level Lindblad model whose transition frequencies
  and dipole strengths come from the Hamiltonian at the actual field and
  strain, solved exactly by full counting statistics, with line
  broadening, spectral diffusion and laser background;
- a 17-check testbench that compares the chain with measurements,
  analytic limits and identities.

## 2. The full file tree

```
snv-tfln-cavity-interface/
|-- README.md            overview, installation, quick start
|-- DETAILS.md           this file
|-- CHANGELOG.md         what changed between versions
|-- CITATION.cff         citation details ("Cite this repository" button)
|-- LICENSE, NOTICE      Apache-2.0 license; copyright and data credits
|-- .zenodo.json         metadata for Zenodo
|-- requirements.txt     Python packages (exact versions)
|-- sim/
|   |-- materials.py         refractive indices from the files in sim/data/
|   |-- waveguide.py         waveguide modes (FEM) and the nominal taper EME
|   |-- run_taper_eme.py     four-mode taper EME of the original run (log for Fig. 2(b))
|   |-- taper_tolerance.py   six-mode taper EME with offset, tilt, gap, sidewalls, slab
|   |-- cavity.py            nanobeam geometry and GME helpers
|   |-- run_cavity.py        one GME cavity solve (Q, V, resonance)
|   |-- run_sweeps.py        GME convergence sweeps (runs run_cavity.py)
|   |-- run_cavity_fields.py production GME solve: field maps, mirror decay
|   |-- run_tolerance.py     GME fabrication tolerances and suspended beam
|   |-- cavity_fdtd.py       3D FDTD of the loaded one-sided cavity (Meep)
|   |-- fdtd_jobs.sh         the FDTD runs used in the article, in order
|   |-- fdtd_analysis.py     Q_L, Purcell factor, power fractions -> cavity_loaded.json
|   |-- BUILD_MEEP.md        how Meep and MPB were built
|   |-- snv.py               SnV- Hamiltonian (levels, cyclicity, qubit frequency, Rabi rate)
|   |-- snv_results.py       spin-model results (validation and scans)
|   |-- readout_fcs.py       four-level Lindblad model and full counting statistics
|   |-- readout_revision.py  all readout predictions of the article -> readout_rev.json
|   |-- checks.py            17-check testbench -> verification.json
|   |-- make_numbers.py      every quoted number as a LaTeX macro -> macros.tex
|   `-- data/                Peter.yml, Zelmon-o.yml, Zelmon-e.yml, Malitson.yml (CC0)
`-- figures/
    |-- figstyle.py          figure sizes, colors, fonts
    |-- fig1_device.py       Fig. 1 (device, level diagram, transition spectrum)
    |-- make_figures.py      Figs. 2 to 5
    `-- fig1_device.png ... fig5_readout.png   rendered previews
```

## Fonts

If you want the figure text in Times New Roman, place the font files
(`times.ttf`, `timesbd.ttf`, `timesi.ttf`, `timesbi.ttf`) in
`figures/fonts/`. They are proprietary and not distributed here. Without
them, matplotlib uses a fallback serif font.

## 4. Files the figure and number scripts expect

All in `sim/results/` (all are in the `data/` folder of the Zenodo
archive):

- waveguide and taper: `wg_diamond.json`, `wg_ln.json`, `wg_hybrid.json`,
  `taper_eme.json`, `eme_log2.txt`, `tolerance_eme/*.json`;
- cavity: `cavity_mode_a*.json`, `cavity_production.json`,
  `cavity_field.npz`, `tolerance.json`, `band_edges.json` (optional; if
  missing, `make_figures.py` recomputes it);
- loaded cavity: `cavity_fdtd_N*_r*.json`, `cavity_fdtd_bulk_r40.json`,
  `cavity_fdtd_summary.json`, `cavity_loaded.json`;
- spin model and readout: `snv_results.json`, `readout_rev.json`;
- testbench: `verification.json` (written by `checks.py`, read by
  `make_numbers.py`).

---

## 5. The scripts, step by step

Run times were measured on one core of a shared two-core computer and are
indicative only.

| Step | Command | What it does | Time | Results (in `sim/results/`) |
|---|---|---|---|---|
| 1 | `python sim/materials.py` | Refractive indices at 619 nm | < 1 s | printed |
| 2 | `python sim/waveguide.py` | Mode studies of the diamond beam, TFLN ridge and hybrid stack | minutes | `wg_diamond.json`, `wg_ln.json`, `wg_hybrid.json` |
| 3 | `python sim/run_taper_eme.py` | Four-mode taper EME of the original submission; its log gives the supermode branches of Fig. 2(b). Needs `gmsh` | about 65 min | `taper_eme.json`, `eme_log2.txt` |
| 4 | `python sim/taper_tolerance.py nominal off50 ...` | Six-mode EME for the nominal taper and each non-ideal geometry (names in `CONFIGS`). Run one configuration per process; the longest sweeps need several GB of memory | 4 to 9 min per configuration | `tolerance_eme/<name>.json` |
| 5 | `python sim/run_cavity.py`, `cd sim; python run_sweeps.py`, `python sim/run_cavity_fields.py`, `python sim/run_tolerance.py` | GME cavity: single solve, convergence sweeps, production run, tolerances | minutes to hours | `cavity_mode_*.json`, `cavity_production.json`, `cavity_field.npz`, `tolerance.json` |
| 6 | `cd sim; bash fdtd_jobs.sh` | 3D FDTD of the loaded cavity for N<sub>out</sub> = 0 to 3 (resolution 40) and N<sub>out</sub> = 0 at resolution 50, plus the bulk reference. Needs Meep | about 25 to 30 min per run at resolution 40 (about 63 min for N<sub>out</sub> = 0, which runs until the field has decayed) and about 1 h at resolution 50 | `cavity_fdtd_*.json` |
| 7 | `python sim/fdtd_analysis.py` | Loaded Q, Purcell factor, power fractions, design point | seconds | `cavity_fdtd_summary.json`, `cavity_loaded.json` |
| 8 | `cd sim; python snv_results.py` | Spin-model validation and scans (must be started inside `sim/`) | seconds | `snv_results.json` |
| 9 | `python sim/readout_revision.py` | All readout predictions (parts `reference`, `ql`, `alpha`, `strain`, `field`, `robust`, `extra`; give part names to run only some) | about 11 min | `readout_rev.json` |
| 10 | `python sim/checks.py` | The 17-check testbench | about 10 s | `verification.json` |
| 11 | `python sim/make_numbers.py` | Writes every number quoted in the article as a LaTeX macro | a few seconds | `macros.tex` |
| 12 | `python figures/fig1_device.py`, `python figures/make_figures.py` | Fig. 1 and Figs. 2 to 5 | about 1 min | `figures/*.pdf`, `figures/*.png` |

Use `OMP_NUM_THREADS=1` (and `OPENBLAS_NUM_THREADS=1`) for the FEM and
readout scripts; with several BLAS threads on a small machine the mode
solves become much slower.

---

## 6. Which script makes which figure

| Figure | Content | Data | Script |
|---|---|---|---|
| Fig. 1 | (a) Top and side views and magnified cavity region (schematic); (b) level diagram with all four C-line transitions and the D line; (c) relative excitation strength of the four C-line transitions versus detuning from A1 | `readout_rev.json` | `fig1_device.py` |
| Fig. 2 | (a) Effective indices of the diamond beam and TFLN ridge; (b) hybrid supermodes along the taper; (c) taper transfer versus length (nominal, 70&deg; sidewalls, ridge on 60-nm slab); (d) lateral offset; (e) angular misalignment; (f) air gap and sidewall angle | steps 2 to 4 | `make_figures.py` |
| Fig. 3 | (a) Band structure and TE gap; (b) cavity field maps; (c) field envelope along the beam; (d) loaded Q versus output mirror periods (FDTD and Fabry-Perot estimate); (e) fate of the emitted power; (f) LDOS spectrum at the design point | steps 5 to 7 | `make_figures.py` |
| Fig. 4 | (a) Cyclicity versus field angle with the measured points; (b) cyclicity versus misalignment; (c) cyclicity and Rabi rate versus strain | step 8 | `make_figures.py` |
| Fig. 5 | (a) Efficiency budget, confocal versus this design; (b) fidelity versus readout window; (c) fidelity versus field misalignment; (d) fidelity versus loaded Q; (e) fidelity and Rabi rate versus strain; (f) fidelity versus line broadening | step 9 | `make_figures.py` |

---

## 7. The Python modules

| File | What it contains |
|---|---|
| `sim/materials.py` | Reads the dispersion formulas from the YAML files in `sim/data/` and evaluates the refractive index; sets the working wavelength (619 nm) |
| `sim/waveguide.py` | femwell/scikit-fem mode solver for the diamond beam, the TFLN ridge and the stack; studies A to C; the four-mode EME of the original submission |
| `sim/taper_tolerance.py` | Six-mode EME on a shared mesh whose diamond layer is divided into 5-nm strips, so that offset, tilt, gap, sidewall angle and a residual slab are described exactly; memory-lean sequential propagation |
| `sim/cavity.py` | legume GME model of the nanobeam: mirror cell, bands, quadratic taper, supercell, Q / V / Purcell factor, TE band-edge classification |
| `sim/cavity_fdtd.py` | Meep model of the open one-sided device: dipole source, Harminv ring-down, flux box, eigenmode decomposition into the waveguide TE<sub>0</sub> mode, LDOS, and the mode volume sampled from the ringing field |
| `sim/fdtd_analysis.py` | Loaded Q, Purcell factor F<sub>cav</sub> = 3Q(&lambda;/n)<sup>3</sup>/(4&pi;<sup>2</sup>V), power fractions at resonance, non-resonant LDOS background and the waveguide fraction &eta;<sub>wg</sub> |
| `sim/snv.py` | SnV<sup>-</sup> Hamiltonian of Rosenthal et al., PRX 13, 031022 (2023), Appendix B, with the parameters of their Table I; dipole operators, cyclicity, qubit frequency, Rabi rate |
| `sim/readout_fcs.py` | Level structure (frequencies and polarization-summed dipole strengths of A1, B2, A2, B1), transition-resolved Purcell budget, Lindblad superoperator, full-counting-statistics generating function, optimal threshold, spectral-diffusion averaging |
| `sim/readout_revision.py` | Reference operating point (9&deg; misalignment, field magnitude fixed by the measured qubit frequency) and every readout scan of the article |
| `figures/figstyle.py` | Figure sizes (full width 379.42 pt, the text width of the optica-article class), colors, fonts |

---

## 8. Where the numbers come from

**Materials** (`sim/data/`, refractiveindex.info, CC0): diamond from
F. Peter, Z. Phys. 15, 358 (1923); LiNbO<sub>3</sub> from Zelmon, Small
and Jundt, J. Opt. Soc. Am. B 14, 3319 (1997); fused silica from
Malitson, J. Opt. Soc. Am. 55, 1205 (1965).

**Waveguide geometry**: diamond beam 200 nm thick (280 nm wide in the
cavity section), following the diamond-on-TFLN platform of Riedel et al.,
ACS Photonics 10, 4236 (2023), whose LN waveguides are ~190-nm ridges on a
~60-nm slab. Our nominal bus is a fully etched 190-nm-high, 600-nm-wide
ridge with vertical sidewalls; the Riedel profile and sloped sidewalls are
evaluated in `taper_tolerance.py`. The taper runs from 350 to 50 nm.

**Cavity design**: hole radius 0.30a, quadratic taper of the hole spacing
over 8 cells to 0.86a (Quan and Loncar, Opt. Express 19, 18529 (2011)),
10 mirror periods on the closed side, mirror period a = 194.35 nm.

**Spin model** (`snv.py`): Hamiltonian and parameters of Rosenthal et
al., PRX 13, 031022 (2023), Appendix B and Table I (spin-orbit 830.15 and
2988.0 GHz, strain 177.67 and 134.00 GHz, g<sub>L</sub> = 0.363 and
0.581, dipole polar angle 125.3&deg; fixed and azimuth 37.33&deg; fitted
in that work), with the ab initio reduction factors of Thiering and Gali,
PRX 8, 021063 (2018). Measured values used for validation: Rosenthal et
al., PRX 14, 041008 (2024) (qubit frequency 3.677 GHz at 125 mT and
&zeta; = 147&deg;; cyclicity 2244 &plusmn; 108 at &zeta; = 147&deg; and
8.6 &plusmn; 0.4 at &zeta; = 53&deg;, for which the field magnitude is not
stated; the model cyclicity is evaluated at 180 mT and changes by less than
1% between 112 and 180 mT).

**Emitter constants** (`readout_fcs.py`): lifetime 4.5 ns and quantum
efficiency 0.80 (Rosenthal et al. 2024); Debye-Waller factor 0.57
(G&ouml;rlitz et al., New J. Phys. 22, 013048 (2020)); C/D branching ratio
0.75 and polarization overlap 2/3 (Lee et al., Phys. Rev. Applied 25,
044074 (2026)).

**Detection chain** (`readout_fcs.cavity_budget`): chip-to-fiber
transmission 0.5, the 3-dB loss of a TFLN edge coupler at 775 nm (Liu et
al., Adv. Photonics Nexus 1, 016001 (2022)); detector efficiency 0.77, the
system detection efficiency of an SNSPD at 610 nm (Wang et al., IEEE
Photonics J. 14, 6816304 (2022)).

**Values set by hand** (device choices and assumptions):

| Value | Used |
|---|---|
| Emitter placement overlap &xi;<sub>pos</sub> | 0.5 (scanned 0.1 to 1) |
| State-preparation fidelity f<sub>0</sub> | 0.99 |
| Detector dark-count rate | 100 per second |
| Field misalignment of the reference point | 9&deg; (reproduces the measured maximum cyclicity) |
| Confocal comparison | &eta; = 0.2% and 0.4%, background 0.2 counts per 50 &micro;s (Rosenthal et al. 2024) |

---

## 9. Built-in checks

`python sim/checks.py` runs 17 checks and writes
`sim/results/verification.json`. Each check prints PASS or FAIL from its
own result, and the last line gives the number passed.

| Check | What it compares | Pass condition |
|---|---|---|
| V1 | Diamond and SiO<sub>2</sub> index at 619 nm against 2.414 and 1.4575 | 0.5%, 0.2% |
| V2 | Zero-field splittings (closed form and diagonalization) against 902.98 and 3000 GHz | 0.2% |
| V3 | Qubit frequency at 125 mT, &zeta; = 147&deg;, against the measured 3.677 GHz | 2% |
| V4 | Cyclicity at &zeta; = 53&deg; (model at 180 mT) against the measured 8.6 | 10% |
| V5 | Cyclicity at 9&deg; misalignment (model at 180 mT) against the measured 2244 | 10% |
| V6 | Spin-flip strength for a field along the spin axis | below 1e-20 |
| V7 | Perturbative qubit frequency (Eq. B9 of Rosenthal 2023) against diagonalization | 1e-3 |
| V8 | Intrinsic Q of the GME convergence runs | largest/smallest < 1.15 |
| V9 | Taper EME: T &le; 1 everywhere; 6 vs 8 modes, 30 vs 45 segments, 5 vs 2.5-nm strips | differences < 0.2% |
| V10 | Four-level model at the confocal point of Rosenthal 2024 | 3 < n<sub>b</sub> < 5, n<sub>d</sub> < 0.3, measured 87.4% between the &eta; = 0.2% and 0.4% predictions |
| V11 | Purcell budget identities: &beta; from transition rates, F<sub>C</sub> - 1 = 4g<sup>2</sup>/(&kappa;&gamma;<sub>C</sub>); bad-cavity condition | exact; g/&kappa; < 0.05 |
| V12 | Lowest intrinsic Q under fabrication perturbations against Q<sub>L</sub> | > 100 Q<sub>L</sub> |
| V13 | Full-counting-statistics normalization | 1e-8 |
| V14 | Two-level limit against the analytic flip-terminated Poisson mixture | CDF within 1% |
| V15 | Mean counts against the time-integrated Lindblad emission | 1e-3 |
| V16 | Readout independent of the field azimuth | 1e-6 |
| V17 | FDTD Purcell factor from the mode volume against the LDOS fit (N<sub>out</sub> = 0, complete ring-down) | 10% |

---

## 10. Notes on the calculations

- **Units.** Waveguide and FDTD geometry in micrometers (Meep time unit
  &micro;m/c). In legume, lengths are in units of the mirror period and
  frequencies in c/a. The spin and readout models use GHz (cycles per ns),
  rates in 1/ns and times in ns.
- **Lithium niobate** is treated with the ordinary index for the TE modes
  of the x-cut film.
- **Taper EME.** Power that leaves the retained modes is lost
  (non-unitary projection), so T &le; 1 (check V9). The overlap is the
  conjugated cross-power integral of femwell. The four-mode EME of the
  original submission (99.1%) under-resolved the TE/TM crossing; the
  six-mode value is used in the revised article. `taper_tolerance.py`
  keeps only two segments in memory and returns freed memory after each
  segment; run one configuration per process.
- **Loaded cavity.** Fractions of the emitted power are evaluated at the
  resonance from flux and mode-monitor spectra of the same ringing field,
  so they do not need a complete ring-down. The Purcell factor uses the
  mode volume sampled from the ringing field; for N<sub>out</sub> = 0 the
  ring-down is complete and the LDOS fit is an independent check (V17).
- **Readout model.** The cavity linewidth (about 0.6 THz) exceeds all
  C-line splittings, so all four C transitions share one Purcell factor and
  the cyclicities are those of the bare emitter. The D line, one
  ground-state splitting away, is enhanced according to the cavity
  Lorentzian. Background is Poissonian; spectral diffusion is averaged by
  nine-point Gauss-Hermite quadrature before the threshold is chosen.
- **Operating point.** The published parameter set does not reproduce both
  the measured qubit frequency and the measured maximum cyclicity at the
  nominal 125 mT, so the reference point uses 9&deg; misalignment and the
  field magnitude (112 mT) at which the model qubit frequency equals
  3.677 GHz. The model also overestimates the A1-B2 splitting at the
  operating point of Rosenthal et al. (2023); the field-magnitude scan
  covers the measured 0.5 GHz.
- **Rabi rate.** `snv.py` prints microwave Rabi rates of about 2 to 3 MHz
  (0.6 mT drive at 125 mT) next to the measured 6.25 MHz of Rosenthal et
  al. (2024); the article uses the model Rabi rate only to show the trend
  with strain.

---

## 11. Version history

| Version | Date | Main change | Archive DOI |
|---|---|---|---|
| v2.0.1 | 2026-10-20 | Figure labels moved clear of the data (Figs. 1, 4, 5) and panel letters aligned; descriptions of checks V4 and V5 state that the source gives no field magnitude for the measured cyclicities; no change in any result | https://doi.org/10.5281/zenodo.23507013 (v2.0.0 code, same results) |
| v2.0.0 | 2026-10-20 | Revised manuscript: four-level readout model, 3D FDTD of the loaded cavity, taper tolerances, measured efficiencies, 17-check testbench; superseded two-level model removed | https://doi.org/10.5281/zenodo.23507013 |
| v1.1.1 | 15 Sep 2026 | Figures resized to the full text width | https://doi.org/10.5281/zenodo.22758086 |
| v1.1.0 | 15 Sep 2026 | New article title; Fig. 1 redrawn; Times New Roman figures | https://doi.org/10.5281/zenodo.22758086 |
| v1.0.0 | 15 Sep 2026 | First version | https://doi.org/10.5281/zenodo.22756033 |
