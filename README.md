# snv-tfln-cavity-interface

Open-source design pipeline for an integrated tin-vacancy (SnV-) spin-photon
interface: an overcoupled diamond nanobeam photonic-crystal cavity,
adiabatically coupled to a thin-film lithium niobate (TFLN) bus, evaluated all
the way to predicted single-shot spin-readout fidelity.

Companion code for the article:

> T. M. Mahim, M. M. Rahman, and A. S. M. Mohsin, "Overcoupled diamond
> nanobeam cavities on thin-film lithium niobate for fast single-shot readout
> of tin-vacancy spins" (submitted to Optics Express, 2026).

The raw simulation outputs, the validation testbench with reference logs, and
a frozen copy of this code are archived on Zenodo:
https://doi.org/10.5281/zenodo.XXXXXXXX

## What the pipeline does

1. **Materials** (`sim/materials.py`) - evaluates CC0-licensed dispersion data
   from refractiveindex.info (diamond, LiNbO3 o/e, SiO2), shipped verbatim in
   `sim/data/`.
2. **Waveguide + taper** (`sim/waveguide.py`) - full-vector FEM cross-section
   modes (femwell / scikit-fem) of the diamond-on-TFLN system at 619 nm, and
   an eigenmode-expansion (EME) computation of the adiabatic taper transfer
   (99.1% at 6 um into a 600-nm LN bus).
3. **Cavity** (`sim/cavity.py`, `sim/run_cavity*.py`) - guided-mode-expansion
   (legume) design of the diamond nanobeam cavity on SiO2:
   Q_i = 1.4e6, V = 0.68 (lambda/n)^3 at 618.7 nm, mirror decay 4.0 dB/period,
   fabrication tolerances, suspended comparison.
4. **Spin model** (`sim/snv.py`, `sim/snv_results.py`) - the strained SnV-
   Hamiltonian of Rosenthal et al. [PRX 13, 031022 (2023)] in QuTiP,
   validated with no free parameters against published measurements
   (level structure, qubit frequency, cyclicity vs field angle).
5. **Interface + readout** (`sim/readout.py`, `sim/interface.py`) - the
   detected-photon budget and exact single-shot counting statistics
   (including dark-state switch-on), design-space maps, and the overcoupled
   design point (Q_L ~ 500, F_C ~ 19, eta ~ 65%, F_r = 98.5% in 0.1 us).
6. **Verification** (`sim/checks.py`) - the 15-check testbench
   (all pass; writes `results/verification.json`).
7. **Figures / numbers** (`figures/make_figures.py`, `sim/make_numbers.py`) -
   regenerate every figure and every number quoted in the article.

## Quick start

```bash
pip install -r requirements.txt
python sim/materials.py          # material sanity check
python sim/snv.py                # spin-model validation numbers
python sim/interface.py          # design point + readout predictions
python sim/checks.py             # full 15-check testbench
python figures/make_figures.py   # all figures (needs results/, see Zenodo)
```

The heavy solves (GME cavity, FEM/EME taper) take minutes to tens of minutes
each; their outputs are archived in the Zenodo dataset so that the figures
and the testbench can be regenerated without re-running them. Place the
archived `data/` contents in `sim/results/` to reuse them.

## Requirements

Python 3.11 with numpy, scipy, matplotlib, qutip, legume-gme, autograd,
femwell, scikit-fem, shapely, gmsh (see `requirements.txt` for pinned
versions). Note: gmsh needs system OpenGL libraries (libglu1-mesa).

## License

Apache-2.0 (see LICENSE, NOTICE). The material data files in `sim/data/` are
CC0 (public domain) from refractiveindex.info.
