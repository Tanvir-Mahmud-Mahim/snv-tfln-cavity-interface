# Fast Single-Shot Readout of Tin-Vacancy Spins: Design and Simulation Code

[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758086.svg)](https://doi.org/10.5281/zenodo.22758086)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

Code for the article **"Fast single-shot readout of tin-vacancy spins with an
overcoupled diamond nanocavity on thin-film lithium niobate"**
by T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin (BRAC University),
submitted to Optics Express (2026).

- Repository: https://github.com/Tanvir-Mahmud-Mahim/snv-tfln-cavity-interface
- Archive on Zenodo (raw simulation outputs, the validation testbench with
  reference logs, and a frozen copy of this code):
  https://doi.org/10.5281/zenodo.22758086

The figure scripts and `sim/make_numbers.py` read most of their numbers
from the saved simulation results in `sim/results/`; the values that are
written directly into these scripts are listed in
[Section 10](#10-notes-on-the-calculations).

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What is in this repository](#2-what-is-in-this-repository)
3. [Installation](#3-installation)
4. [Quick start: three ways to use the code](#4-quick-start-three-ways-to-use-the-code)
5. [The scripts, step by step](#5-the-scripts-step-by-step)
6. [Which script makes which figure](#6-which-script-makes-which-figure)
7. [The Python modules](#7-the-python-modules)
8. [Where the numbers come from](#8-where-the-numbers-come-from)
9. [Built-in checks](#9-built-in-checks)
10. [Notes on the calculations](#10-notes-on-the-calculations)
11. [Version history](#11-version-history)
12. [How to cite](#12-how-to-cite)
13. [License and contact](#13-license-and-contact)

---

## 1. The idea in one minute

A **tin-vacancy centre** (written SnV<sup>-</sup>) is a tin atom trapped in a
gap of missing carbon atoms in diamond. It carries an electron spin that can
store a quantum bit, and it gives off red light near 619 nm.

To read the spin, a laser makes the centre glow when the spin points one way
(the "bright" state) and stay nearly dark when it points the other way. If
enough photons reach the detector, one counting window is enough to tell the
two apart: this is **single-shot readout**. Two things limit it. Few photons
reach the detector, and every emitted photon carries a small chance of
flipping the spin, which ends the signal. The average number of photons a
centre emits before its spin flips is called the **cyclicity**.

This project designs a chip that collects the light better and faster:

- a narrow diamond bar (a **nanobeam**) with a row of holes that forms an
  optical **cavity**, a tiny resonator that traps light next to the centre;
- the cavity makes the centre emit faster and mostly into the cavity (the
  **Purcell effect**; the speed-up is the **Purcell factor**);
- the cavity is **overcoupled**: it is built to leak its light into the
  output waveguide much faster than it loses it anywhere else;
- a gradually narrowing end of the diamond bar (a **taper**) hands the light
  to a waveguide in **thin-film lithium niobate** (TFLN), which carries it
  across the chip.

The code follows the light from start to finish:

- refractive indices of diamond, lithium niobate and glass (SiO<sub>2</sub>);
- the light modes of the diamond bar and the TFLN waveguide, and how much
  light the taper transfers from one to the other;
- the cavity: its colour, how long it stores light (the quality factor
  **Q**) and how tightly it squeezes the light (the mode volume **V**);
- a model of the SnV<sup>-</sup> spin checked against published measurements
  (energy levels, qubit frequency, cyclicity versus magnet angle);
- the photon budget and the probability of reading the spin correctly in one
  shot (the **readout fidelity** F<sub>r</sub>), with a 15-item testbench
  that checks the whole chain.

---

## 2. What is in this repository

```
snv-tfln-cavity-interface/
|-- README.md               this guide
|-- CHANGELOG.md            what changed between versions
|-- CITATION.cff            citation details (drives the "Cite this repository" button)
|-- LICENSE                 Apache-2.0 license
|-- NOTICE                  copyright notice and credit for the refractiveindex.info data
|-- .zenodo.json            metadata for Zenodo (title, authors, keywords)
|-- requirements.txt        Python packages to install (exact versions)
|-- sim/
|   |-- materials.py        refractive indices read from the files in sim/data/
|   |-- waveguide.py        light modes of the waveguides (FEM) and the taper transfer (EME)
|   |-- cavity.py           nanobeam cavity geometry and solver helpers (GME)
|   |-- run_cavity.py       one cavity solve: finds the cavity mode, Q, V and Purcell factor
|   |-- run_sweeps.py       runs run_cavity.py for the mirror-number and convergence sweeps
|   |-- run_cavity_fields.py  final ("production") cavity solve: field maps and mirror decay
|   |-- run_tolerance.py    fabrication-tolerance and suspended-beam cavity solves
|   |-- snv.py              SnV- spin model (levels, qubit frequency, cyclicity, Rabi rate)
|   |-- snv_results.py      spin-model results for the article (validation and scans)
|   |-- readout.py          photon budget and single-shot readout statistics
|   |-- interface.py        design point, readout predictions and the design-space map
|   |-- checks.py           the verification testbench (writes verification.json)
|   |-- make_numbers.py     collects every quoted number into LaTeX macros
|   `-- data/               four refractiveindex.info files (CC0): Peter.yml (diamond),
|                           Zelmon-o.yml, Zelmon-e.yml (LiNbO3), Malitson.yml (SiO2)
`-- figures/
    |-- figstyle.py         shared figure style (sizes, Okabe-Ito colours, fonts)
    |-- fig1_device.py      draws Fig. 1 (device schematic and headline result)
    |-- make_figures.py     draws Figs. 2 to 5
    `-- fig1_device.png ... fig5_readout.png   rendered previews of the five figures
```

The simulation scripts write their results to `sim/results/`
(`snv_results.py` writes to `results/` in the folder it is started from).
The folder is created automatically and is not stored on GitHub (it is
listed in `.gitignore`); the archived results are in the Zenodo record. The figure scripts write PDF and
PNG files into `figures/`. `sim/make_numbers.py` writes `paper/macros.tex`
and, if the files `paper/manuscript.tex` and `supplement/supplement.tex`
exist (they are not in this repository), also replaces the block between
their "auto-generated result macros" markers with the new macros.

---

## 3. Installation

You need **Python 3.11 or newer**: the pinned `numpy`, `scipy` and `qutip`
versions do not install on older Python. The previous README named
Python 3.11, and the checks below were run with Python 3.11.15.

```
pip install -r requirements.txt
```

This installs fixed versions of `numpy`, `scipy`, `matplotlib`, `qutip` (the
spin model), `legume-gme` and `autograd` (the cavity solver), `femwell`,
`scikit-fem`, `shapely`, `gmsh` and `meshio` (the waveguide solver). `pip`
also installs what these need, for example `meshwell` and `pygmsh` for
`femwell`, in versions that are not pinned.

**System library.** `gmsh` needs the OpenGL library `libGLU.so.1` (package
`libglu1-mesa` on Debian and Ubuntu). Without it, `import gmsh` fails and
`sim/waveguide.py` cannot run. All other scripts work without it.

**Fonts (optional).** Figure text uses Times New Roman when its font files
(`times.ttf`, `timesbd.ttf`, `timesi.ttf`, `timesbi.ttf`) are placed in
`figures/fonts/`. They are not distributed here; without them matplotlib
uses a fallback font.

---

## 4. Quick start: three ways to use the code

Run all commands from the repository folder unless a step says otherwise.

### Way A: run the fast models (no downloads, under a minute)

```
python sim/materials.py      # refractive indices at 619 nm
python sim/snv.py            # spin model next to the measured values
python sim/interface.py      # design point and readout fidelity
```

`sim/interface.py` prints the design point (loaded Q<sub>L</sub> = 500,
Purcell factor F<sub>C</sub> of about 19.3, total detection efficiency of
about 65%) and the best readout fidelity for three values of the bare
cyclicity. For the highest cyclicity (2244) it prints F<sub>r</sub> = 0.985
with a 0.08 &micro;s window (quoted as "98.5% in 0.1 &micro;s" in the
earlier README and in Fig. 1).

`python sim/checks.py` also runs without downloads, but only its first ten
checks (V1 to V7); it then stops with an error, because checks V8 to V12
read the saved cavity, taper and tolerance results (Way B).

### Way B: full testbench, figures and numbers from the archived results

1. Download the archive from https://doi.org/10.5281/zenodo.22758086.
2. Put the contents of its `data/` folder into `sim/results/` (this is the
   instruction given in the earlier README; the file list of the Zenodo
   record could not be opened while this guide was written).
3. Run:

```
python sim/checks.py             # full testbench -> sim/results/verification.json
python figures/fig1_device.py    # Fig. 1
python figures/make_figures.py   # Figs. 2 to 5
python sim/make_numbers.py       # LaTeX macros -> paper/macros.tex
```

The scripts expect these files in `sim/results/`: `wg_diamond.json`,
`wg_ln.json`, `taper_eme.json`, `eme_log2.txt`, `cavity_mode_a193*.json`,
`cavity_production.json`, `cavity_field.npz`, `tolerance.json`,
`snv_results.json`, `interface.json`, `design_map.npz` and
`fr_vs_lambda.json` (plus `verification.json` for `make_numbers.py`, which
`checks.py` writes). `band_edges.json` is optional: if it is missing,
`make_figures.py` computes it.

### Way C: recompute the results yourself

Run the steps in [Section 5](#5-the-scripts-step-by-step) in order. The
cavity and waveguide solves are the slow part; the earlier README gives
"minutes to tens of minutes" for each heavy solve. Three result files are
**not** produced by any command in this repository:

- `taper_eme.json` and `eme_log2.txt` come from the function
  `eme_taper()` in `sim/waveguide.py`, which no script calls. The arguments
  used for the archived run (bus width, taper lengths) are not recorded.
  The default length list is 2, 4, 6, 8, 10, 14 and 18 &micro;m, but
  `make_numbers.py` needs an entry at 3 &micro;m, and `make_figures.py`
  labels the fifth entry as the 6 &micro;m value, so the archived run used a
  different list. `make_figures.py` reads `eme_log2.txt` as the printed log
  of such a run (lines of the form `[i/30] w=... n_eff=[...]`).
- `fr_vs_lambda.json` (readout fidelity versus bare cyclicity, for the cavity
  design and for a confocal set-up with detection efficiency 0.2% and 0.4%)
  is read by `fig1_device.py`, `make_figures.py` and `make_numbers.py`, but
  no script here writes it.

For these files, use the archived copies (Way B).

---

## 5. The scripts, step by step

| Step | Command | What it does | Time* | Results (in `sim/results/`) |
|---|---|---|---|---|
| 1 | `python sim/materials.py` | Prints the refractive indices of diamond, LiNbO<sub>3</sub> (ordinary and extraordinary) and SiO<sub>2</sub> at 619 nm, and the diamond group index | 0.4 s | printed only |
| 2 | `python sim/waveguide.py` | Light modes (effective index versus width) of a 200 nm thick diamond bar on SiO<sub>2</sub> (widths 160 to 400 nm) and of a 190 nm thick TFLN ridge (widths 250 to 700 nm) | not run here (needs `libGLU`) | `wg_diamond.json`, `wg_ln.json` |
| 2b | *(no script)* `eme_taper()` in `sim/waveguide.py` | Light transfer through the diamond taper (350 nm down to 50 nm, 30 slices) onto the TFLN bus, versus taper length | not run here | `taper_eme.json`; the printed log is `eme_log2.txt` |
| 3 | `python sim/run_cavity.py [a gmax numeig depth N_mirror W_y]` | One cavity solve: finds the mode trapped at the centre, then its Q, mode volume and Purcell factor. Defaults: a = 180 nm, gmax 2.0, 130 modes, taper depth 0.14, 10 mirror holes, W_y = 4 | 3.7 min (223 s) with the defaults | `cavity_mode_a<a>_g<gmax>_d<depth>_N<N>_W<W_y>.json` |
| 4 | `cd sim` then `python run_sweeps.py` | Runs `run_cavity.py` at a = 193.4 nm for 4, 6, 8 and 12 mirror holes, for gmax 1.5, 2.5 and 3.0, and for W_y = 5 (8 solves). Must be started inside `sim/`, because it calls `run_cavity.py` by its bare file name | long (not re-timed) | `cavity_mode_a193_*.json` |
| 5 | `python sim/run_cavity_fields.py` | Final cavity solve (a = 194.35 nm, gmax 2.5): wavelength, Q, field maps, and the light decay per mirror period | long (not re-timed) | `cavity_production.json`, `cavity_field.npz` |
| 6 | `python sim/run_tolerance.py` | Seven cavity solves at a = 194.35 nm: nominal, hole radius 0.29a and 0.31a, width 270 and 290 nm, thickness 190 nm, and a free-standing (suspended) beam | long (not re-timed) | `tolerance.json` |
| 7 | `cd sim` then `python snv_results.py` | Spin-model numbers: checks against measurements, cyclicity versus magnet angle and misalignment, and cyclicity and Rabi rate versus strain. Must be started inside `sim/`: it writes to `results/` relative to the current folder | 6 s | `snv_results.json` |
| 8 | `python sim/interface.py --map` | Design point (Q<sub>L</sub> = 500), readout for three bare cyclicities, best Q<sub>L</sub> scans, and the fidelity map over Q<sub>L</sub> and emitter placement. Without `--map` the map is skipped (33 s) | 106 s | `interface.json`, `design_map.npz` |
| 9 | `python sim/checks.py` | The verification testbench ([Section 9](#9-built-in-checks)) | 3 s for V1 to V7 | `verification.json` |
| 10 | `python figures/fig1_device.py` and `python figures/make_figures.py` | Draw Fig. 1 and Figs. 2 to 5 | not timed (needs archived data) | `figures/*.pdf`, `figures/*.png` |
| 11 | `python sim/make_numbers.py` | Writes every number quoted in the article as a LaTeX macro; also splices them into `paper/manuscript.tex` and `supplement/supplement.tex` if those files exist | not timed (needs archived data) | `paper/macros.tex` |

\*Times measured on a shared two-core computer that was busy with other work
(load average about 6).
`python sim/readout.py` (0.5 s) prints a check of the confocal readout model
and a demonstration case; it saves nothing.

`run_sweeps.py` does not include a run at gmax 2.0 with 10 mirror holes and
W_y = 4 at a = 193.4 nm; `python sim/run_cavity.py 193.4` makes that one.
Check V8 and `make_numbers.py` use every `cavity_mode_a193*.json` file with
10 mirror holes that is present.

Steps 2 to 6 (waveguide and cavity solves) are slow. Steps 1 and 7 to 11
are fast.

---

## 6. Which script makes which figure

`make_figures.py` loads all of its input files before drawing, so every file
listed in Way B must be present even to draw one figure. It always draws
Figs. 2 to 5 in one run.

| Figure in the article | Content | Data from | Drawn by |
|---|---|---|---|
| Fig. 1 | (a) Top view, side view and magnified cavity region of the device (schematic); (b) level scheme with the cavity-enhanced transition; (c) readout fidelity versus bare cyclicity | `fr_vs_lambda.json` | `fig1_device.py` |
| Fig. 2 | (a) Effective index versus width, diamond bar and TFLN ridge; (b) coupled-mode effective indices along the taper; (c) taper transfer versus taper length | steps 2, 2b | `make_figures.py` (`fig2`) |
| Fig. 3 | (a) Band diagram of one mirror cell with the TE band gap and the cavity line (computed while drawing); (b) field maps of the cavity mode, top and side; (c) field along the beam on a log scale | step 5 (and `band_edges.json`) | `make_figures.py` (`fig3`) |
| Fig. 4 | (a) Cyclicity versus magnet angle, with the two measured points; (b) cyclicity versus misalignment; (c) cyclicity and Rabi rate versus strain | step 7 | `make_figures.py` (`fig4`) |
| Fig. 5 | (a) Efficiency of each stage, confocal set-up versus this design; (b) readout fidelity versus bare cyclicity; (c) fidelity map over Q<sub>L</sub> and emitter placement, with the validity boundary | step 8 and `fr_vs_lambda.json` | `make_figures.py` (`fig5`) |

---

## 7. The Python modules

| File | What it contains |
|---|---|
| `sim/materials.py` | Reads the dispersion formulas and coefficients from the YAML files in `sim/data/` (no retyped numbers) and evaluates the refractive index; refuses wavelengths outside each file's range. Sets the working wavelength, 619 nm |
| `sim/waveguide.py` | Cross-section mode solver (femwell, finite elements on gmsh meshes) for the diamond bar, the TFLN ridge and the two stacked; studies A to C; eigenmode-expansion (EME) transfer through the taper |
| `sim/cavity.py` | Guided-mode-expansion (legume) model of the nanobeam: mirror cell, band structure, quadratic taper of the hole spacing, cavity supercell, Q / mode volume / Purcell factor, and `classify_te_edges()` for the TE band edges |
| `sim/snv.py` | Spin Hamiltonian of the strained SnV<sup>-</sup> ground and excited states (spin-orbit, strain, Zeeman and orbital Zeeman terms), dipole operators, cyclicity, qubit frequency and microwave Rabi rate |
| `sim/readout.py` | Emission budget constants, Purcell and detuning formulas, readout photon counts, and single-shot fidelity models (`fidelity`, `fidelity_threshold`, `fidelity_exact`) including the dark-state "switch-on" process |
| `sim/interface.py` | Links the cavity, taper and spin results: `design()` gives all derived quantities at one loaded Q; `readout_point()` finds the best fidelity over threshold, window length and laser strength |
| `figures/figstyle.py` | Figure sizes (single column 85 mm; full width 379.42 pt, the text width of the optica-article LaTeX class), colours and fonts |

---

## 8. Where the numbers come from

The values below are written in the code, with the sources the code gives.

**Materials** (`sim/data/`, copied from refractiveindex.info, CC0):
diamond from F. Peter, Z. Phys. 15, 358 (1923); LiNbO<sub>3</sub> (ordinary
and extraordinary) from Zelmon, Small and Jundt, J. Opt. Soc. Am. B 14, 3319
(1997); fused silica from Malitson, J. Opt. Soc. Am. 55, 1205 (1965). At
619 nm `materials.py` gives n = 2.4137 (diamond), 2.2903 (LiNbO<sub>3</sub>,
ordinary), 2.2055 (extraordinary) and 1.4574 (SiO<sub>2</sub>). The working
wavelength 619 nm is the SnV<sup>-</sup> zero-phonon line (the code cites
Rugar et al., Phys. Rev. B 99, 205417 (2019), and Rosenthal et al., PRX 14,
041008 (2024)).

**Waveguide geometry** (`waveguide.py`) follows the diamond-on-TFLN platform
of Riedel et al. (arXiv:2306.15207, doi:10.1364/OPTICAQ.1.000103): diamond
200 nm thick, lithium-niobate film 190 nm thick on SiO<sub>2</sub>. The
chosen widths are 280 nm (diamond bar) and 600 nm (TFLN bus), as written in
`make_numbers.py` and Fig. 1.

**Cavity design** (`cavity.py`, `run_cavity_fields.py`): hole radius
r = 0.30a; hole spacing reduced quadratically over 8 cells on each side of
the centre to 0.86a (taper depth 0.14), following Quan and Loncar, Opt.
Express 19, 18529 (2011); 10 mirror holes on each side; final mirror period
a = 194.35 nm, basis size gmax 2.5.

**Spin model** (`snv.py`): Hamiltonian of Rosenthal et al., Phys. Rev. X
13, 031022 (2023), Appendix B, with values from its Table I: spin-orbit
splittings 830.15 GHz (ground) and 2988.0 GHz (excited); strain 177.67 GHz
and 134.00 GHz; factors g<sub>L</sub> = 0.363 and 0.581, with reduction
factors from Thiering and Gali, Phys. Rev. X 8, 021063 (2018). The dipole
direction uses a polar angle of 125.3&deg;, marked "fixed" in the `snv.py`
docstring, and an azimuthal angle of 37.33&deg;, marked "(fit)"; the
docstring does not say where this fit was made. The electron gyromagnetic
ratio is 28 GHz/T. The measured values used for
comparison (qubit frequency 3.677 GHz at 125 mT; cyclicity 2244 and 8.6 at
180 mT) are cited in `snv.py` as Rosenthal et al., arXiv:2403.13110; check
V3 in `checks.py` cites the qubit frequency as PRX 14, 041008.

**Emitter constants** (`readout.py`): excited-state lifetime 4.5 ns and
quantum efficiency 0.80 (arXiv:2403.13110), Debye-Waller factor 0.57
(Goerlitz et al. 2020; Lee et al., arXiv:2511.05740), C/D branching ratio
0.75 (Lee et al., arXiv:2511.05740), dipole alignment factor sin<sup>2</sup>(54.7&deg;) = 2/3.

**Photonic inputs to the readout model** (constants written in
`interface.py`; its comments say they come from the GME cavity and EME taper
results): cavity wavelength 618.9 nm, intrinsic Q 1.39e6,
mode volume 0.68 (&lambda;/n)<sup>3</sup>, mirror intensity decay per period
0.399, taper transfer 0.991 at 6 &micro;m, qubit frequency 3.677 GHz, C/D
splitting 2.097 THz, laser saturation parameter 2.

*Conflict, not resolved here:* the previous README gave the cavity as
"Q_i = 1.4e6, V = 0.68 (lambda/n)^3 at 618.7 nm, mirror decay 4.0
dB/period", while `interface.py` uses 618.9 nm (its comment reads
"production GME; tuned to C line"). The mirror values agree: a decay factor
of 0.399 per period equals 3.99 dB per period, which rounds to 4.0 dB. The
production wavelength itself is saved in `cavity_production.json` (Zenodo
archive), which could not be opened to check which value it holds.

**Values set in the code without a cited source** (device choices and
assumptions):

| Value | Used |
|---|---|
| Loaded quality factor Q<sub>L</sub> (design point) | 500 |
| Emitter placement overlap &xi;<sub>pos</sub> | 0.5 (the map scans 0.05 to 1.0) |
| On-chip routing and off-chip coupling efficiency | 0.8 |
| Detector efficiency | 0.9 |
| State-preparation fidelity f<sub>0</sub> | 0.99 |
| Background count rate | 1000 per second |
| Mirror loading model: group index, cavity length | 10, 1.2 &micro;m |

---

## 9. Built-in checks

`sim/checks.py` runs 12 groups of checks and records 15 results in
`sim/results/verification.json`, ending with a line "n/15 checks passed".

| Check | What it compares | Allowed deviation |
|---|---|---|
| V1a, V1b | Diamond and SiO<sub>2</sub> index at 619 nm against 2.414 and 1.4575 | 0.5%, 0.2% |
| V2a to V2c | Zero-field splittings (formula and full calculation) against 902.98 GHz and 3000 GHz | 0.2% |
| V3 | Qubit frequency at 125 mT, 147&deg; against the measured 3.677 GHz | 2% |
| V4 | Cyclicity at 53&deg;, 180 mT against the measured 8.6 | 10% |
| V5 | Cyclicity with the field 9&deg; off the spin axis against the measured maximum 2244 | 10% |
| V6 | With the field exactly along the spin axis, the spin-flip strength must vanish (below 1e-20) | - |
| V7 | Approximate qubit-frequency formula (Eq. B9) against full calculation at 50 mT | 2% |
| V8 | Intrinsic Q from all saved `cavity_mode_a193*.json` runs with 10 mirror holes (different gmax and W_y) | largest/smallest below 1.15 |
| V9 | Taper transfer never above 1, and the last two lengths differ by less than 0.001 | - |
| V10 | Readout model with confocal detection efficiency 0.2 to 0.4% against the measured case (about 4 bright counts; measured fidelity 87.4%) | about 3 to 5 counts; fidelity not above 88% |
| V11 | Purcell factor from Q and V equals that from the coupling rates (an identity) | 1e-6 |
| V12 | Lowest intrinsic Q in the tolerance runs must exceed 100 &times; 500; the wavelength shifts are recorded | - |

For V8, V9 and V12 the script always prints the word "PASS"; the real
pass/fail result is the `ok` field saved in `verification.json`, which is
what the final count uses.

`sim/snv.py` and `sim/readout.py` also print their results next to the
measured values when run on their own.

---

## 10. Notes on the calculations

- **Units.** Waveguide geometry is in micrometres. Inside the cavity solver
  (legume) lengths are in units of the mirror period a and frequencies in
  units of c/a; `mode_metrics()` converts to nanometres and
  (&lambda;/n)<sup>3</sup>. The spin model works in GHz (energies divided by
  Planck's constant) and tesla.
- **Lithium niobate** is treated as isotropic with the ordinary index in
  the waveguide solves (the code comment says the ordinary and
  extraordinary values are quoted in the manuscript).
- **Taper transfer (EME).** Light that leaves the kept set of modes is
  counted as lost (the code calls this a non-unitary projection); check V9
  tests that the transfer never exceeds 100%. All slices share one mesh so
  that the mode overlaps are computed exactly.
- **Cavity mode choice.** The run scripts pick, among the modes in the band
  gap, the one with the largest share of its light within 3 periods of the
  centre. `run_cavity_fields.py` uses a fixed band-edge value (0.2884 c/a);
  its comment says this is the value at a = 193.4 nm from a `run_cavity.py`
  output, while the script itself solves a = 194.35 nm.
- **Cyclicity at the measured maximum.** With the magnet angle of the
  measurement (147&deg;, 180 mT) the model gives a cyclicity of about 158
  (`snv.py` prints 158.5 next to the measured 2244). The testbench (V5) and
  the strain scan instead use a field 9&deg; away from the spin axis, which
  gives about 2193; the code note refers to a reported residual
  misalignment of about 10&deg;.
- **Rabi rate.** `snv.py` prints microwave Rabi rates of 2.75, 3.27 and
  1.89 MHz (drive angles 0, 45 and 90&deg;; 0.6 mT drive at 125 mT) next to
  the measured 6.25 MHz, so the model gives roughly a third to a half of the
  measured value. No check in `checks.py` compares the Rabi rate.
- **Readout statistics.** `interface.py` uses `fidelity_exact()`: bright
  counts are a Poisson process that stops when the spin flips; dark counts
  include background, weak off-resonant light, and the chance that the dark
  spin flips into the bright state during the window. The best integer
  photon threshold is chosen for each case.
- **Two different example cases.** The demonstration at the bottom of
  `readout.py` uses Q = 1e5, V = 0.73 and 85% out-coupling; the design point
  of the article is the one in `interface.py` (Q<sub>L</sub> = 500). The
  docstring at the top of `interface.py` still mentions an older design
  point (Q<sub>L</sub> about 1.5e3, F<sub>C</sub> about 57); the code runs
  Q<sub>L</sub> = 500.
- **Values written directly into the output scripts.** Fig. 5(a) takes its
  per-stage efficiencies from numbers written in `make_figures.py`. Among
  them is the cavity fraction &beta; = 0.905, while `interface.py`
  computes &beta; = 0.909 (0.90925) at the design point; the other entries
  for this design (0.9996, 0.991, 0.8, 0.9) match the code. Fig. 1
  writes "a = 194 nm", "Q<sub>L</sub> &asymp; 500", "F<sub>C</sub> = 19" and
  "98.5% in 0.1 &micro;s" as text; `make_numbers.py` writes the mode volume
  (0.68, spread 0.65 to 0.78), the widths (280 nm, 600 nm), the taper length
  (6 &micro;m), Q<sub>L</sub> = 500 and the confocal efficiency range
  (0.2 to 0.4%) directly.
- **Band diagram of Fig. 3(a)** is recomputed with legume each time
  `make_figures.py` runs.
- **Output message.** `interface.py` prints "saved interface.json +
  design_map.npz" even without `--map`, when only `interface.json` is
  written.

---

## 11. Version history

Dates are commit dates of the tagged versions.

| Version | Date | Main change | Archive DOI cited in the README |
|---|---|---|---|
| v1.1.1 (latest tag) | 15 Sep 2026 | Figures resized to the full text width (379.4 pt) so text prints at its nominal size | https://doi.org/10.5281/zenodo.22758086 |
| v1.1.0 | 15 Sep 2026 | New article title; Fig. 1 redrawn as 2D views (`fig1_device.py`); Times New Roman figures; TE band-edge classification added | https://doi.org/10.5281/zenodo.22758086 |
| v1.0.0 | 15 Sep 2026 | First version: simulation code, material data, figure scripts | https://doi.org/10.5281/zenodo.22756033 |

The `main` branch has four more commits after v1.1.1 (figure label fixes,
15 to 16 Sep 2026), and this documentation follows them; neither is tagged.
`CITATION.cff` therefore still gives version 1.1.1 dated 2026-09-15, the
last tag. Both are listed under "Unreleased" in [CHANGELOG.md](CHANGELOG.md).

---

## 12. How to cite

Please cite the article and the Zenodo archive. GitHub also shows a
**"Cite this repository"** button in the right-hand column, which reads
`CITATION.cff`.

> T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin, "Fast single-shot readout
> of tin-vacancy spins with an overcoupled diamond nanocavity on thin-film
> lithium niobate" (submitted to Optics Express, 2026).
>
> Archive: https://doi.org/10.5281/zenodo.22758086

---

## 13. License and contact

Code: Apache License 2.0 (see `LICENSE` and `NOTICE`). The material data
files in `sim/data/` are from the refractiveindex.info database and are in
the public domain (CC0).

Questions and bug reports: please open an issue on this repository, or
contact Tanvir M. Mahim, BRAC University (tanvir.mahim@bracu.ac.bd).
