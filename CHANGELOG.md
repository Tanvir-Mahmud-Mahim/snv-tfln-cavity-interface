# Changelog

All notable changes to this code are listed here, newest first. Entries
before "Unreleased" are taken from the git history and tags.

## Unreleased

### Documentation

Documentation only; no code, data or result changed.

- README rewritten as a step-by-step guide: the idea in plain language, a
  file tree, installation (Python 3.11 or newer; the `libGLU` system library
  needed by `gmsh`), three ways to use the code, a table of the scripts with
  run times measured on a shared two-core computer, which script makes which
  figure, where the numbers come from, what the testbench compares, and
  notes on the calculations.
- The README now states which result files no script in the repository
  writes (`taper_eme.json`, `eme_log2.txt`, `fr_vs_lambda.json`) and that
  `run_sweeps.py` and `snv_results.py` must be started inside `sim/`.
- Added this CHANGELOG.
- CITATION.cff: added the version, release date and Zenodo DOI; the article
  status is now given as `status: submitted`. The version is 1.1.1 dated
  2026-09-15, the last tag; the figure commits below and this documentation
  come after that tag and are not part of any tagged version.

### Figures (15 to 16 September 2026, after v1.1.1)

- Figure label clashes fixed after checks at high zoom (Figs. 1 to 5),
  in commits `d52c7c6`, `73237af`, `d64b66b` and `97fd8f6`. The y-axis of
  Fig. 2(a) now extends below the SiO2 light line, and the two-line note in
  Fig. 4(a) became a short "see (b)" label (the commit message says the
  explanation now lives in the manuscript caption).

## v1.1.1 (15 September 2026)

- README wording fix: the Fig. 1 schematic is 2D, not 3D (`c8a6000`).
- Figures drawn at the full text width of the optica-article LaTeX class
  (379.4 pt, `figures/figstyle.py`), so that figure text prints at its
  nominal point size. All figure PNGs regenerated.

## v1.1.0 (15 September 2026)

- New article title: "Fast single-shot readout of tin-vacancy spins with an
  overcoupled diamond nanocavity on thin-film lithium niobate" (README and
  CITATION.cff).
- Zenodo DOI in the README changed to
  https://doi.org/10.5281/zenodo.22758086 (the commit message calls it the
  dataset DOI for v1.1).
- Fig. 1 redesigned as 2D top, side and magnified views, drawn by the new
  script `figures/fig1_device.py`; the old `figures/fig1_pipeline.png` was
  removed. Figure label clashes fixed in Figs. 2 to 5.
- Figure text set in Times New Roman when the font files are placed in
  `figures/fonts/`.
- Added `classify_te_edges()` to `sim/cavity.py` (TE band-edge
  classification used by Fig. 3).
- Corrected the zero-phonon-line reference comment in `sim/materials.py`.

## v1.0.0 (15 September 2026)

First version.

- Repository root: README, Apache-2.0 LICENSE, NOTICE, CITATION.cff,
  `.zenodo.json`, `requirements.txt`.
- `sim/`: materials, FEM/EME waveguide, GME cavity, spin model, readout
  model, design integration and the verification testbench.
- `sim/data/`: CC0 dispersion files from refractiveindex.info (diamond,
  LiNbO3 ordinary and extraordinary, SiO2).
- `figures/`: figure scripts and rendered PNG previews.
- Article title at this version: "Overcoupled diamond nanobeam cavities on
  thin-film lithium niobate for fast single-shot readout of tin-vacancy
  spins". The README cited the reserved Zenodo DOI
  https://doi.org/10.5281/zenodo.22756033.
