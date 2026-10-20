# Changelog

All notable changes to this code are listed here, newest first. Entries
for v1.x are taken from the git history and tags.

## v2.0.1 (2026-10-20)

Presentation fixes only; every number, figure curve and check result is
unchanged.

- `figures/make_figures.py`: Fig. 5(b) labels moved clear of the curves;
  the omega_e = omega_q label of Fig. 5(c) sits above the frame; the Fig.
  5(f) legend no longer covers the last data point; panel letters in a row
  share one height; Fig. 4 uses the same letter placement.
- `figures/fig1_device.py`: the loaded-Q label of the magnified view sits
  under the output section, inside the panel; the A1 label of panel (c)
  is moved clear of its curve.
- `sim/checks.py`, `sim/snv.py`: the descriptions of checks V4 and V5 no
  longer say that the measured cyclicities were taken at 180 mT; the source
  does not state the field magnitude, and the model value (computed at
  180 mT) changes by less than 1% between 112 and 180 mT.

## v2.0.0 (2026-10-20)

Code for the revised manuscript (revision 1, Optics Express). The
archived data are in https://doi.org/10.5281/zenodo.23507013.

### New

- `sim/readout_fcs.py`: four-level Lindblad model of the readout
  (states |dn>, |up>, |A>, |B>). The laser drives all four C-line
  transitions (A1, B2, A2, B1) with Rabi frequencies scaled by the square
  root of each transition's polarization-summed dipole strength and
  detuned by each transition's frequency offset, all computed from the
  Hamiltonian in `snv.py`. The detected-count distributions are computed
  exactly by full counting statistics (generating function on the unit
  circle, inverted by FFT). Extra homogeneous broadening (pure dephasing),
  quasi-static spectral diffusion (Gauss-Hermite average before the
  threshold is chosen) and Poissonian laser background are included.
- `sim/readout_revision.py`: every readout prediction of the revised
  article (reference operating point, fidelity versus window, loaded Q,
  field misalignment, strain, field magnitude, line broadening,
  background, detection chain, placement overlap, D-line detection,
  direct fiber coupling, azimuthal invariance) -> `readout_rev.json`.
- `sim/taper_tolerance.py`: six-mode taper EME on a shared mesh whose
  diamond layer is divided into 5-nm strips, for lateral offset, angular
  misalignment, a diamond-LN air gap, sloped LN sidewalls and a partially
  etched ridge on a 60-nm slab -> `tolerance_eme/*.json`.
- `sim/cavity_fdtd.py`, `sim/fdtd_jobs.sh`, `sim/fdtd_analysis.py`,
  `sim/BUILD_MEEP.md`: 3D FDTD (Meep) of the loaded one-sided cavity:
  loaded Q (Harminv), Purcell factor from the mode volume of the ringing
  field (checked against the LDOS), and the fractions of the emission
  coupled into the output waveguide TE0 mode, the substrate, air, the
  sides and the closed mirror -> `cavity_fdtd_summary.json`,
  `cavity_loaded.json`.

### Changed

- `sim/checks.py`: 17 checks (V1 to V17). New checks of the readout model
  (V10 confocal operating point, V13 normalization, V14 two-level limit,
  V15 mean counts against the Lindblad emission, V16 azimuthal
  invariance), the taper convergence (V9: 6 vs 8 modes, 30 vs 45 segments,
  strip width) and the FDTD Purcell factor (V17). V3 and V5 now cite the
  field of each measurement. Every check prints PASS or FAIL from its own
  result (before, V8, V9 and V12 always printed PASS).
- `sim/make_numbers.py`: rewritten for the revised article; every number
  comes from a results file (no hand-typed design values); writes
  `sim/results/macros.tex`.
- `figures/fig1_device.py`: panel (b) is the complete level diagram (all
  four C-line transitions and the D line); panel (c) shows the relative
  excitation strength of the four C-line transitions versus detuning from
  A1; the loaded Q label is read from the FDTD result.
- `figures/make_figures.py`: Fig. 2 has six panels (taper transfer for
  realistic LN profiles, lateral offset, tilt, gap and sidewall angle);
  Fig. 3 adds the FDTD loaded Q, the fate of the emitted power and the
  LDOS spectrum; Fig. 5 is rebuilt from `readout_rev.json` (efficiency
  budget, fidelity versus window, misalignment, loaded Q, strain and line
  broadening).
- `sim/waveguide.py`: the docstring cited a wrong DOI for Riedel et al.
  and called their LN layer a 190-nm film; it is a ~190-nm ridge on a
  ~60-nm slab (ACS Photonics 10, 4236 (2023)).
- `sim/snv.py`: the docstring now says where the dipole angles come from
  (Table I of Rosenthal et al., PRX 13, 031022 (2023): polar angle fixed,
  azimuth fitted).
- `requirements.txt`: adds h5py (needed by Meep) and a note on Meep.

### Removed

- `sim/readout.py`, `sim/interface.py`, `sim/run_fr_vs_lambda.py`: the
  two-level readout model of the original submission. It admitted
  dark-state excitation only through the spin-flipping transition and
  with the strength of the cycling transition, omitted the dark spin's
  spin-conserving transition B2, placed the D line 2.097 THz from the C
  line (it is one ground-state splitting, about 903 GHz, away), used
  beta = eta_q eta_DW eta_BR F_C / zeta instead of (zeta - 1)/zeta, and
  assumed chip and detector efficiencies of 0.8 and 0.9. These files
  remain in release v1.1.1.

### Results that changed

- Taper transfer at 6 um: 99.7% with six retained modes (the four-mode
  value of v1, 99.1%, under-resolved the TE/TM crossing).
- Loaded cavity: Q_L = 870 (FDTD) with 92% of the cavity decay in
  the output waveguide (v1: Q_L = 500 assumed, 99.96% estimated).
- Detected-photon efficiency 32% (v1: 65%); single-shot fidelity
  98.8% in 0.1 us (v1: 98.5%).

## Changes after v1.1.1 (untagged until v2.0.0)

### Fixes (30 September 2026)

No scientific result changed. The regenerated files were compared with the
archived ones in the Zenodo record (10.5281/zenodo.22758086).

- New script `sim/run_taper_eme.py`. Before: `taper_eme.json` and
  `eme_log2.txt` came from `eme_taper()` in `sim/waveguide.py`, which no
  script called, and the settings of the archived run were not recorded.
  After: the script calls `eme_taper()` with the settings recovered from
  the archived files (bus width 0.6 um, taper 350 nm to 50 nm in 30 slices,
  4 modes per slice, lengths 1, 2, 3, 4, 6, 8, 10, 12, 16 and 20 um) and
  saves the printed log as `eme_log2.txt`. Result: the log is byte-for-byte
  identical to the archived `eme_log2.txt`; the transfer values in
  `taper_eme.json` differ from the archived ones by at most 3.6e-14.
- New script `sim/run_fr_vs_lambda.py`. Before: `fr_vs_lambda.json` (read by
  `fig1_device.py`, `make_figures.py` and `make_numbers.py`) was written by
  no script. After: the script writes it; the cavity curve uses
  `interface.readout_point()` at Q_L = 500, and the confocal curves use
  `readout.fidelity_exact()` with the check-V10 settings, state-preparation
  fidelity 0.99 and the better of a 10 us and a 30 us window. These
  confocal settings were recovered by comparison with the archived file.
  Result: byte-for-byte identical to the archived `fr_vs_lambda.json`.
- `sim/make_numbers.py`. Before: it created `paper/` in the repository
  folder, wrote `paper/macros.tex`, and replaced the macro block inside
  `paper/manuscript.tex` and `supplement/supplement.tex` when those files
  existed. After: it writes the same macros to `sim/results/macros.tex`
  (the folder is created if needed and is already in `.gitignore`) and no
  longer creates folders or edits any manuscript file. The macro text is
  unchanged (checked byte for byte against the old script).
- README: the two new scripts in the file tree, in Way C and in the script
  table (steps 2b and 8b), how they were checked, the new place of
  `macros.tex`, and the measured run times of `make_numbers.py` and the
  figure scripts. The "no script writes these files" gap is removed.

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
