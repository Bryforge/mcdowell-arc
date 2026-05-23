import numpy as np

from mcdowell_arc.constants import OMEGA_EARTH_RAD_S, R_EARTH_KM
from mcdowell_arc.frames import earth_rotation_velocity_km_s, geodetic_to_ecef_km
from mcdowell_arc.orbit import summarize_orbit


def test_geodetic_equator_radius():
    r = geodetic_to_ecef_km(0, 0, 0)
    assert np.allclose(r, [R_EARTH_KM, 0, 0])


def test_rotation_velocity_equator():
    r = np.array([R_EARTH_KM, 0, 0])
    v = earth_rotation_velocity_km_s(r)
    assert np.allclose(v, [0, OMEGA_EARTH_RAD_S * R_EARTH_KM, 0])


def test_circular_orbit_perigee_close_to_altitude():
    alt = 400.0
    r = R_EARTH_KM + alt
    v = np.sqrt(398600.4418 / r)
    summary = summarize_orbit(np.array([r, 0, 0, 0, v, 0]))
    assert abs(summary.perigee_km - alt) < 1e-6
    assert abs(summary.apogee_km - alt) < 1e-6
