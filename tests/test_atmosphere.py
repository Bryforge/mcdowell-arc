import numpy as np

from mcdowell_arc.atmosphere import density_kg_m3


def test_density_monotonically_decreases_in_upper_atmosphere():
    samples = density_kg_m3(np.array([100.0, 150.0, 250.0, 400.0]))
    assert np.all(np.diff(samples) < 0.0)


def test_density_scalar_return_type():
    assert isinstance(density_kg_m3(100.0), float)
