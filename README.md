# Fast Single-Shot Readout of Tin-Vacancy Spins: Design and Simulation Code

[![Data DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23507013.svg)](https://doi.org/10.5281/zenodo.23507013)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

This repository holds the code for our article **"Fast single-shot readout
of tin-vacancy spins with an overcoupled diamond nanocavity on thin-film
lithium niobate"**, which I wrote with M. M. Rahman and A. S. M. Mohsin
(BRAC University). It has been submitted to Optics Express (2026); this
version (v2.0.1) belongs to the revised manuscript (revision 1). The Zenodo
archive holds v2.0.0, which differs only in figure label positions and
check descriptions and gives identical results.

- Repository: https://github.com/Tanvir-Mahmud-Mahim/snv-tfln-cavity-interface
- Archive on Zenodo (raw simulation outputs, the validation testbench with
  reference outputs, and a frozen copy of this code):
  https://doi.org/10.5281/zenodo.23507013

Every number quoted in the revised article is written by
`sim/make_numbers.py` from the saved simulation results in
`sim/results/`. The few values that are set by hand (device choices and
measured literature values) are listed in
[Section 8 of DETAILS.md](DETAILS.md#8-where-the-numbers-come-from).

---

## Contents

1. [The idea in one minute](#1-the-idea-in-one-minute)
2. [What changed in version 2.0.0](#2-what-changed-in-version-200)
3. [What is in this repository](#3-what-is-in-this-repository)
4. [Installation](#4-installation)
5. [Quick start](#5-quick-start)
6. [How to cite](#6-how-to-cite)
7. [License and contact](#7-license-and-contact)

The step-by-step guide to every script, the figures, the modules, the
sources of all numbers and the testbench are in [DETAILS.md](DETAILS.md).

---

## 1. The idea in one minute

A **tin-vacancy center** (SnV<sup>-</sup>) is a tin atom trapped in a gap
of missing carbon atoms in diamond. It carries an electron spin that can
store a quantum bit, and it gives off red light near 619 nm.

To read the spin, you shine a laser on the center. The center glows when
the spin points one way (the "bright" state) and stays nearly dark when it
points the other way. If enough photons reach the detector, one counting
window tells the two apart: this is **single-shot readout**. Two things
limit it: few photons reach the detector, and every emitted photon carries
a small chance of flipping the spin. A third effect matters as soon as the
emission is collected efficiently: the laser also excites the dark spin
through its own spin-conserving transition, a few hundred megahertz to a
few gigahertz away, and those photons look like bright-state photons.

We design a chip that collects the light efficiently without paying that
price:

- a diamond **nanobeam** with a row of holes forms an optical **cavity**
  that makes the center emit faster and mostly into the cavity (the
  **Purcell effect**);
- the cavity is **overcoupled**: one of its mirrors is shortened so that it
  leaks into the output waveguide quickly. Its moderate Purcell factor keeps
  the broadened optical line narrow compared with the splitting between the
  bright and dark transitions;
- a **taper** at the end of the diamond beam hands the light to a waveguide
  in **thin-film lithium niobate** (TFLN), which routes it off the chip.

The main results of the revised article, as `sim/make_numbers.py` writes
them:

- loaded cavity (3D FDTD): Q<sub>L</sub> = 870, Purcell factor
  F<sub>C</sub> = 28, 92% of the cavity decay reaches the output
  waveguide;
- taper transfer 99.7% at 6 &micro;m, and 99.4% for a combined
  realistic placement and etch error;
- detected-photon efficiency 32% with measured chip-to-fiber (0.5)
  and detector (0.77) efficiencies;
- single-shot readout fidelity 98.8% in a 0.1 &micro;s window from a
  four-level model computed from the SnV<sup>-</sup> Hamiltonian at the
  actual field and strain, above 98.1% for every strain studied.

---

## 2. What changed in version 2.0.0

Version 2.0.0 accompanies the revised manuscript. It adds the simulations
that the reviewers asked for and corrects errors that we found while
revising. The full list is in [CHANGELOG.md](CHANGELOG.md). In short:

- **New four-level readout model** (`sim/readout_fcs.py`,
  `sim/readout_revision.py`): a Lindblad master equation for the two
  ground and two excited spin states, solved exactly by full counting
  statistics, with every transition frequency and dipole strength taken
  from the Hamiltonian. It includes the dark spin's own transition B2,
  dephasing, spectral diffusion and laser background.
- **Loaded cavity in 3D FDTD** (`sim/cavity_fdtd.py`,
  `sim/fdtd_analysis.py`, Meep): Q<sub>L</sub>, Purcell factor and where
  the emitted power goes (waveguide, substrate, air, closed mirror).
- **Taper tolerances** (`sim/taper_tolerance.py`): lateral offset,
  angular misalignment, diamond-LN gap, sloped sidewalls and a partially
  etched ridge.
- **Testbench** (`sim/checks.py`): now 17 checks (V1 to V17), each printed
  as PASS or FAIL from its own result.
- **Corrections**: the D line is placed one ground-state splitting
  (about 903 GHz) from the C line, not 2.097 THz; the cavity fraction is
  &beta; = 1 - &gamma;<sub>0</sub>/&gamma;<sub>cav</sub>; the taper EME
  keeps six modes instead of four; the chip-to-fiber and detector
  efficiencies are measured values instead of assumptions.
- **Removed**: the two-level readout model of the original submission
  (`sim/readout.py`, `sim/interface.py`, `sim/run_fr_vs_lambda.py`). It is
  superseded and contained the errors listed above. It remains available
  in release [v1.1.1](https://github.com/Tanvir-Mahmud-Mahim/snv-tfln-cavity-interface/releases/tag/v1.1.1)
  and in the earlier Zenodo versions.

---

## 3. What is in this repository

- `sim/`: simulation modules and scripts; `sim/data/` holds four
  refractiveindex.info dispersion files (CC0).
- `sim/BUILD_MEEP.md`: how we built Meep and MPB for the FDTD runs.
- `figures/`: figure scripts and rendered previews of Figs. 1 to 5.
- `requirements.txt`: Python packages (exact versions).
- [DETAILS.md](DETAILS.md): the full guide (file tree, every script, run
  times, figures, sources of numbers, testbench).
- [CHANGELOG.md](CHANGELOG.md): what changed between versions.
- `CITATION.cff`, `.zenodo.json`: citation metadata.
- `LICENSE`, `NOTICE`: Apache-2.0 license and credits.

The scripts write their results to `sim/results/`. That folder is not on
GitHub (it is in `.gitignore`); the archived results are in the Zenodo
record. The manuscript and the supplement are not part of this repository.

---

## 4. Installation

You need **Python 3.11 or newer**. I ran the revision with Python 3.13.

```
pip install -r requirements.txt
```

`gmsh` needs the system library `libGLU.so.1` (package `libglu1-mesa` on
Debian and Ubuntu); without it the waveguide and taper scripts cannot run.

The FDTD scripts (`cavity_fdtd.py`) need **Meep** with **MPB**, which are
not on PyPI. Install them from conda-forge (`conda install -c conda-forge
pymeep`) or build them from source as described in
[sim/BUILD_MEEP.md](sim/BUILD_MEEP.md). Everything else runs without
Meep, including the readout model, the figures and the testbench, as long
as the archived FDTD results are in `sim/results/`.

How to get Times New Roman in the figures (optional) is in
[DETAILS.md](DETAILS.md#fonts).

---

## 5. Quick start

**Fast models only (no downloads):**

```
python sim/materials.py      # refractive indices at 619 nm
python sim/snv.py            # spin model next to the measured values
```

**Figures, numbers and testbench from the archived results:**

1. Download the archive from https://doi.org/10.5281/zenodo.23507013.
2. Copy the contents of its `data/` folder into `sim/results/`.
3. Run:

```
python sim/checks.py              # 17-check testbench -> verification.json
python sim/make_numbers.py        # LaTeX macros -> sim/results/macros.tex
python figures/fig1_device.py     # Fig. 1
python figures/make_figures.py    # Figs. 2 to 5
```

**Recompute everything:** follow
[Section 5 of DETAILS.md](DETAILS.md#5-the-scripts-step-by-step). The
waveguide, cavity and FDTD solves take minutes to hours each.

---

## 6. How to cite

Please cite the article and the Zenodo archive. GitHub also shows a
**"Cite this repository"** button, which reads `CITATION.cff`.

> T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin, "Fast single-shot readout
> of tin-vacancy spins with an overcoupled diamond nanocavity on thin-film
> lithium niobate" (submitted to Optics Express, 2026).
>
> Archive: https://doi.org/10.5281/zenodo.23507013

---

## 7. License and contact

Code: Apache License 2.0 (see `LICENSE` and `NOTICE`). The material data
files in `sim/data/` are from the refractiveindex.info database and are in
the public domain (CC0).

If you have questions or find a bug, please open an issue on this
repository. You can also contact me, Tanvir M. Mahim, at BRAC University
(tanvir.mahim@bracu.ac.bd).
