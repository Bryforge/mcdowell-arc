"""Frame helpers for a lightweight ECI/ECEF MVP.

This module intentionally avoids high-precision Earth orientation machinery. For
serious analysis, replace the simple rotation with astropy/ERFA/IERS-backed time
handling.
"""

from __future__ import annotations

import numpy as np

from .constants import OMEGA_EARTH_RAD_S, R_EARTH_KM


def rot_z(theta_rad: float) -> np.ndarray:
    """Return a right-handed rotation matrix about +Z."""
    c = np.cos(theta_rad)
    s = np.sin(theta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def geodetic_to_ecef_km(lat_deg: float, lon_deg: float, alt_km: float) -> np.ndarray:
    """Spherical-Earth geodetic coordinates to ECEF position in kilometers."""
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    r = R_EARTH_KM + alt_km
    return np.array(
        [
            r * np.cos(lat) * np.cos(lon),
            r * np.cos(lat) * np.sin(lon),
            r * np.sin(lat),
        ]
    )


def ecef_to_eci_position_km(r_ecef_km: np.ndarray, t_s: float, theta0_rad: float = 0.0) -> np.ndarray:
    """Convert an ECEF position vector to the simplified ECI frame."""
    theta = theta0_rad + OMEGA_EARTH_RAD_S * t_s
    return rot_z(theta) @ np.asarray(r_ecef_km, dtype=float)


def earth_rotation_velocity_km_s(r_eci_km: np.ndarray) -> np.ndarray:
    """Velocity of the co-rotating atmosphere at an ECI position."""
    omega = np.array([0.0, 0.0, OMEGA_EARTH_RAD_S])
    return np.cross(omega, np.asarray(r_eci_km, dtype=float))


def earth_relative_to_inertial_velocity_km_s(
    r_eci_km: np.ndarray, v_rel_eci_km_s: np.ndarray
) -> np.ndarray:
    """Convert an Earth-relative velocity vector to inertial velocity.

    This requires the velocity direction in the same local/inertial basis. If a
    webcast gives only speed magnitude, this conversion is underconstrained.
    """
    return np.asarray(v_rel_eci_km_s, dtype=float) + earth_rotation_velocity_km_s(r_eci_km)
