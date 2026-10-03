# Fast Single-Shot Readout of Tin-Vacancy Spins: Design and Simulation Code

[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22758086.svg)](https://doi.org/10.5281/zenodo.22758086)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

This repository holds the code for our article **"Fast single-shot readout
of tin-vacancy spins with an overcoupled diamond nanocavity on thin-film
lithium niobate"**, which I wrote with M. M. Rahman and A. S. M. Mohsin
(BRAC University). It has been submitted to Optics Express (2026).

- Repository: https://github.com/Tanvir-Mahmud-Mahim/snv-tfln-cavity-interface
- Archive on Zenodo (raw simulation outputs, the validation testbench with
  reference logs, and a frozen copy of this code):
  https://doi.org/10.5281/zenodo.22758086

The figure scripts and `sim/make_numbers.py` read most of their numbers
from the saved simulation results in `sim/results/`. Some values are
written directly into these scripts, and I list them in
[Section 10](DETAILS.md#10-notes-on-the-calculations).

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What is in this repository](#2-what-is-in-this-repository)
3. [Installation](#3-installation)
4. [Quick start: three ways to use the code](#4-quick-start-three-ways-to-use-the-code)
5. [The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step)
6. [Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure)
7. [The Python modules](DETAILS.md#7-the-python-modules)
8. [Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from)
9. [Built-in checks](DETAILS.md#9-built-in-checks)
10. [Notes on the calculations](DETAILS.md#10-notes-on-the-calculations)
11. [Version history](DETAILS.md#11-version-history)
12. [How to cite](#12-how-to-cite)
13. [License and contact](#13-license-and-contact)

Sections 5 to 11 are in [DETAILS.md](DETAILS.md). So are the
[extra notes for Sections 1 to 4](DETAILS.md#extra-notes-for-sections-1-to-4).

---

## 1. The idea in one minute

A **tin-vacancy center** (written SnV<sup>-</sup>) is a tin atom trapped in a
gap of missing carbon atoms in diamond. It carries an electron spin that can
store a quantum bit, and it gives off red light near 619 nm.

To read the spin, you shine a laser on the center. The center glows when
the spin points one way (the "bright" state). It stays nearly dark when
the spin points the other way. If enough photons reach the detector, one
counting window is enough to tell the two apart. This is **single-shot
readout**. Two things limit it. First, few photons reach the detector.
Second, every emitted photon carries a small chance of flipping the spin,
which ends the signal. The average number of photons a center emits
before its spin flips is called the **cyclicity**.

In this project we design a chip that collects the light better and faster:

- a narrow diamond bar (a **nanobeam**) with a row of holes that forms an
  optical **cavity**, a tiny resonator that traps light next to the center;
- the cavity makes the center emit faster and mostly into the cavity (the
  **Purcell effect**; the speed-up is the **Purcell factor**);
- the cavity is **overcoupled**: it is built to leak its light into the
  output waveguide much faster than it loses it anywhere else;
- a gradually narrowing end of the diamond bar (a **taper**) hands the light
  to a waveguide in **thin-film lithium niobate** (TFLN), which carries it
  across the chip.

The code follows the light from start to finish. It covers the
materials, the waveguides and the taper, the cavity, and a model of the
SnV<sup>-</sup> spin. It ends with the photon budget and the readout
fidelity F<sub>r</sub>, with a 15-item testbench that checks the whole
chain. The full list is in
[DETAILS.md](DETAILS.md#more-on-section-1-what-the-code-computes).

The main results, as `sim/interface.py` prints them:

- the design point: loaded Q<sub>L</sub> = 500, Purcell factor
  F<sub>C</sub> of about 19.3, total detection efficiency of about 65%;
- for the highest bare cyclicity (2244), F<sub>r</sub> = 0.985 with a
  0.08 &micro;s window. The earlier README and Fig. 1 quote this as
  "98.5% in 0.1 &micro;s".

---

## 2. What is in this repository

The main parts are:

- `sim/`: the simulation scripts and modules. Its folder `sim/data/` holds
  four refractiveindex.info files (CC0).
- `figures/`: the figure scripts and rendered previews of the five figures.
- requirements.txt: Python packages to install (exact versions).
- [CHANGELOG.md](CHANGELOG.md): what changed between versions.
- `CITATION.cff`: citation details (drives the "Cite this repository" button).
- `LICENSE` and `NOTICE`: the Apache-2.0 license; copyright notice and
  credit for the refractiveindex.info data.

The simulation scripts write their results to `sim/results/`
(`snv_results.py` writes to `results/` in the folder it is started from).
The folder is created automatically. It is not stored on GitHub (it is
listed in `.gitignore`), but you can find the archived results in the
Zenodo record. The figure scripts write PDF and PNG files into
`figures/`.

The full annotated file tree and more notes on the output folders are in
[DETAILS.md](DETAILS.md#more-on-section-2-the-full-file-tree).

---

## 3. Installation

You need **Python 3.11 or newer**, because the pinned `numpy`, `scipy` and
`qutip` versions do not install on older Python. The previous README named
Python 3.11, and I ran my checks with Python 3.11.15.

```
pip install -r requirements.txt
```

This installs fixed versions of:

- `numpy`, `scipy` and `matplotlib`;
- `qutip` (the spin model);
- `legume-gme` and `autograd` (the cavity solver);
- `femwell`, `scikit-fem`, `shapely`, `gmsh` and `meshio` (the waveguide
  solver).

`pip` also installs what these need, for example `meshwell` and `pygmsh`
for `femwell`. Those versions are not pinned.

**System library.** `gmsh` needs the OpenGL library `libGLU.so.1` (package
`libglu1-mesa` on Debian and Ubuntu). Without it, `import gmsh` fails, and
`sim/waveguide.py` and `sim/run_taper_eme.py` cannot run. All other
scripts work without it.

How to get Times New Roman in the figures (optional) is in
[DETAILS.md](DETAILS.md#more-on-section-3-optional-fonts).

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
about 65%). It also prints the best readout fidelity for three values of
the bare cyclicity. For the highest cyclicity (2244) it prints
F<sub>r</sub> = 0.985 with a 0.08 &micro;s window. The earlier README and
Fig. 1 quote this as "98.5% in 0.1 &micro;s".

You can also run `python sim/checks.py` without downloads, but only its
first ten checks (V1 to V7) will run. It then stops with an error,
because checks V8 to V12 read the saved cavity, taper and tolerance
results (Way B).

### Way B: full testbench, figures and numbers from the archived results

1. Download the archive from https://doi.org/10.5281/zenodo.22758086.
2. Put the contents of its `data/` folder into `sim/results/`. This is the
   instruction given in the earlier README. I could not open the file list
   of the Zenodo record while writing this guide.
3. Run:

```
python sim/checks.py             # full testbench -> sim/results/verification.json
python figures/fig1_device.py    # Fig. 1
python figures/make_figures.py   # Figs. 2 to 5
python sim/make_numbers.py       # LaTeX macros -> sim/results/macros.tex
```

The list of files the scripts expect in `sim/results/` is in
[DETAILS.md](DETAILS.md#way-b-the-files-the-scripts-expect).

### Way C: recompute the results yourself

Run the steps in [Section 5](DETAILS.md#5-the-scripts-step-by-step) in order. The
cavity and waveguide solves are the slow part. The earlier README gives
"minutes to tens of minutes" for each heavy solve.

Two small scripts make the three result files that earlier had no script
(`taper_eme.json`, `eme_log2.txt` and `fr_vs_lambda.json`):

```
python sim/run_taper_eme.py      # taper_eme.json and eme_log2.txt (needs gmsh)
python sim/run_fr_vs_lambda.py   # fr_vs_lambda.json
```

How the settings of these two scripts were recovered, and how the results
compare with the archived ones, is in
[DETAILS.md](DETAILS.md#way-c-how-the-settings-of-the-two-small-scripts-were-recovered).

---

## More details (Sections 5 to 11)

The full notes are in [DETAILS.md](DETAILS.md). Here is what each section holds:

- [5. The scripts, step by step](DETAILS.md#5-the-scripts-step-by-step):
  a table of every step, with its command, what it does, its run time and
  its results.
- [6. Which script makes which figure](DETAILS.md#6-which-script-makes-which-figure):
  what each figure shows, its data and the script that draws it.
- [7. The Python modules](DETAILS.md#7-the-python-modules):
  what each module contains.
- [8. Where the numbers come from](DETAILS.md#8-where-the-numbers-come-from):
  the sources of the material data, geometry, cavity design, spin model and
  emitter constants, and the values set without a cited source.
- [9. Built-in checks](DETAILS.md#9-built-in-checks):
  the 15 results of the testbench and what each one compares.
- [10. Notes on the calculations](DETAILS.md#10-notes-on-the-calculations):
  units, model choices, and values written directly into the output scripts.
- [11. Version history](DETAILS.md#11-version-history):
  the tagged versions, their dates and the archive DOIs.

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

If you have questions or find a bug, please open an issue on this
repository. You can also contact me, Tanvir M. Mahim, at BRAC University
(tanvir.mahim@bracu.ac.bd).
