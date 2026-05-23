"""Keplerian summary quantities from state vectors."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import MU_EARTH_KM3_S2, R_EARTH_KM


@dataclass(frozen=True)
class OrbitSummary:
    perigee_km: float
    apogee_km: float | None
    eccentricity: float
    semi_major_axis_km: float | None
    specific_energy_km2_s2: float


def summarize_orbit(state: np.ndarray) -> OrbitSummary:
    """Return conic summary values for an inertial state vector.

    Perigee is computed for elliptic, parabolic, or hyperbolic conics. Apogee is
    only meaningful for bound elliptic cases.
    """
    state = np.asarray(state, dtype=float)
    r = state[:3]
    v = state[3:]
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)
    h_vec = np.cross(r, v)
    h = np.linalg.norm(h_vec)
    if r_norm <= 0.0 or h <= 0.0:
        raise ValueError("Degenerate orbit state.")

    energy = 0.5 * v_norm**2 - MU_EARTH_KM3_S2 / r_norm
    e_vec = np.cross(v, h_vec) / MU_EARTH_KM3_S2 - r / r_norm
    ecc = float(np.linalg.norm(e_vec))
    p = h**2 / MU_EARTH_KM3_S2
    rp = p / (1.0 + ecc)
    perigee_alt = rp - R_EARTH_KM

    if energy < 0.0:
        a = -MU_EARTH_KM3_S2 / (2.0 * energy)
        ra = a * (1.0 + ecc)
        apogee_alt = ra - R_EARTH_KM
    else:
        a = None
        apogee_alt = None

    return OrbitSummary(
        perigee_km=float(perigee_alt),
        apogee_km=None if apogee_alt is None or not math.isfinite(apogee_alt) else float(apogee_alt),
        eccentricity=float(ecc),
        semi_major_axis_km=None if a is None or not math.isfinite(a) else float(a),
        specific_energy_km2_s2=float(energy),
    )
