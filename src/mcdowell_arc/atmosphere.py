"""Simple atmospheric density model for the MVP.

The density table is a rough high-altitude approximation, log-interpolated by
altitude. It is good enough for software plumbing and sensitivity tests; serious
trajectory claims should use a vetted model such as NRLMSISE-00 with solar and
geomagnetic inputs.
"""

from __future__ import annotations

import numpy as np

# altitude_km, density_kg_m3 rough table from sea level through exosphere-ish altitudes.
_DENSITY_TABLE = np.array(
    [
        [0.0, 1.225],
        [25.0, 3.899e-2],
        [30.0, 1.774e-2],
        [40.0, 3.972e-3],
        [50.0, 1.057e-3],
        [60.0, 3.206e-4],
        [70.0, 8.770e-5],
        [80.0, 1.905e-5],
        [90.0, 3.396e-6],
        [100.0, 5.297e-7],
        [110.0, 9.661e-8],
        [120.0, 2.438e-8],
        [130.0, 8.484e-9],
        [140.0, 3.845e-9],
        [150.0, 2.070e-9],
        [180.0, 5.464e-10],
        [200.0, 2.789e-10],
        [250.0, 7.248e-11],
        [300.0, 2.418e-11],
        [350.0, 9.518e-12],
        [400.0, 3.725e-12],
        [450.0, 1.585e-12],
        [500.0, 6.967e-13],
        [600.0, 1.454e-13],
        [700.0, 3.614e-14],
        [800.0, 1.170e-14],
        [900.0, 5.245e-15],
        [1000.0, 3.019e-15],
    ],
    dtype=float,
)


def density_kg_m3(altitude_km: float | np.ndarray) -> float | np.ndarray:
    """Return log-interpolated atmospheric density in kg/m^3."""
    alt = np.asarray(altitude_km, dtype=float)
    alt_clipped = np.clip(alt, _DENSITY_TABLE[0, 0], _DENSITY_TABLE[-1, 0])
    log_rho = np.interp(
        alt_clipped,
        _DENSITY_TABLE[:, 0],
        np.log(_DENSITY_TABLE[:, 1]),
    )
    rho = np.exp(log_rho)
    if np.isscalar(altitude_km):
        return float(rho)
    return rho
