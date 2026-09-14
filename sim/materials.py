"""Refractive-index models from the refractiveindex.info database (CC0).

Each function evaluates the dispersion formula recorded in the CC0 YAML files
stored in ``sim/data/`` (copied verbatim from the refractiveindex.info
database, https://refractiveindex.info, public domain / CC0):

* ``Peter.yml``      -- diamond (C), formula 1 (Sellmeier-2 form),
                        F. Peter, Z. Phys. 15, 358-368 (1923),
                        doi:10.1007/BF01330487; range 0.226-0.760 um.
* ``Zelmon-o.yml``   -- LiNbO3, ordinary ray, formula 2,
                        D. E. Zelmon, D. L. Small, D. Jundt,
                        J. Opt. Soc. Am. B 14, 3319-3322 (1997),
                        doi:10.1364/JOSAB.14.003319; range 0.40-5.0 um.
* ``Zelmon-e.yml``   -- LiNbO3, extraordinary ray, formula 2, same reference.
* ``Malitson.yml``   -- fused silica (SiO2), formula 1,
                        I. H. Malitson, J. Opt. Soc. Am. 55, 1205-1209 (1965),
                        doi:10.1364/JOSA.55.001205; range 0.21-6.7 um.

The coefficients are parsed from the YAML files at import time so that the
numbers used in the simulations are exactly the numbers in the shipped
database files (no retyping).
"""

from __future__ import annotations

import pathlib
import re

import numpy as np

DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

# SnV zero-phonon-line wavelength (um).  The C and D transitions lie
# near 619-620 nm (Rugar et al., Phys. Rev. B 99, 205417 (2019));
# Rosenthal et al., PRX 14, 041008 (2024) operate at 619.140 nm.
LAMBDA_SNV_UM = 0.619


def _load_formula(path: pathlib.Path) -> tuple[int, np.ndarray, tuple[float, float]]:
    """Parse formula number, coefficients and validity range from a
    refractiveindex.info YAML file (formula 1 or 2 records only)."""
    text = path.read_text()
    m = re.search(r"type:\s*formula\s*(\d)", text)
    if not m:
        raise ValueError(f"no formula record in {path}")
    ftype = int(m.group(1))
    m = re.search(r"coefficients:\s*([-\d.eE+\s]+)", text)
    coeffs = np.array([float(x) for x in m.group(1).split()])
    m = re.search(r"wavelength_range:\s*([\d.]+)\s+([\d.]+)", text)
    if m:
        rng = (float(m.group(1)), float(m.group(2)))
    else:
        rng = (0.0, np.inf)
    return ftype, coeffs, rng


def _n_formula(ftype: int, c: np.ndarray, lam_um) -> np.ndarray:
    """Evaluate refractiveindex.info dispersion formula 1 or 2.

    Formula 1 (Sellmeier):     n^2 - 1 = c0 + sum_i c_{2i+1} lam^2/(lam^2 - c_{2i+2}^2)
    Formula 2 (Sellmeier-2):   n^2 - 1 = c0 + sum_i c_{2i+1} lam^2/(lam^2 - c_{2i+2})
    """
    lam2 = np.asarray(lam_um, dtype=float) ** 2
    n2m1 = c[0] * np.ones_like(lam2)
    for i in range(1, len(c) - 1, 2):
        B, C = c[i], c[i + 1]
        denom = lam2 - (C ** 2 if ftype == 1 else C)
        n2m1 = n2m1 + B * lam2 / denom
    return np.sqrt(1.0 + n2m1)


class Material:
    def __init__(self, name: str, yml: str):
        self.name = name
        self.path = DATA_DIR / yml
        self.ftype, self.coeffs, self.range_um = _load_formula(self.path)

    def n(self, lam_um=LAMBDA_SNV_UM):
        lam = np.asarray(lam_um, dtype=float)
        lo, hi = self.range_um
        if np.any(lam < lo) or np.any(lam > hi):
            raise ValueError(
                f"{self.name}: wavelength {lam} um outside tabulated "
                f"range [{lo}, {hi}] um"
            )
        return _n_formula(self.ftype, self.coeffs, lam)

    def group_index(self, lam_um=LAMBDA_SNV_UM, dlam=1e-4):
        lam = float(lam_um)
        dn = (self.n(lam + dlam) - self.n(lam - dlam)) / (2 * dlam)
        return float(self.n(lam)) - lam * dn


diamond = Material("diamond", "Peter.yml")
linbo3_o = Material("LiNbO3 (o)", "Zelmon-o.yml")
linbo3_e = Material("LiNbO3 (e)", "Zelmon-e.yml")
silica = Material("SiO2", "Malitson.yml")


def summary(lam_um=LAMBDA_SNV_UM):
    rows = []
    for mat in (diamond, linbo3_o, linbo3_e, silica):
        rows.append((mat.name, float(mat.n(lam_um))))
    return rows


if __name__ == "__main__":
    lam = LAMBDA_SNV_UM
    print(f"Refractive indices at lambda = {lam*1e3:.1f} nm (SnV ZPL):")
    for name, n in summary(lam):
        print(f"  {name:12s}  n = {n:.4f}")
    print(f"  diamond group index n_g = {diamond.group_index(lam):.4f}")
