import numpy as np
import pytest

from mcdowell_arc.constants import MU_EARTH_KM3_S2, R_EARTH_KM
from mcdowell_arc.dynamics import propagate
from mcdowell_arc.fast import core_available


def test_propagate_python_returns_expected_shape():
    r = R_EARTH_KM + 400.0
    v = np.sqrt(MU_EARTH_KM3_S2 / r)
    states = propagate(np.array([0.0, 10.0, 20.0]), np.array([r, 0, 0, 0, v, 0]), 250.0, use_drag=False, backend="python")
    assert states.shape == (3, 6)
    assert np.all(np.isfinite(states))


def test_auto_backend_falls_back_without_rust():
    r = R_EARTH_KM + 400.0
    v = np.sqrt(MU_EARTH_KM3_S2 / r)
    states = propagate(np.array([0.0, 5.0]), np.array([r, 0, 0, 0, v, 0]), 250.0, use_drag=False, backend="auto")
    assert states.shape == (2, 6)


@pytest.mark.skipif(not core_available(), reason="Rust extension is not installed")
def test_rust_backend_matches_python_for_short_no_drag_propagation():
    r = R_EARTH_KM + 400.0
    v = np.sqrt(MU_EARTH_KM3_S2 / r)
    times = np.array([0.0, 10.0])
    state0 = np.array([r, 0, 0, 0, v, 0])
    py_states = propagate(times, state0, 250.0, use_drag=False, backend="python")
    rs_states = propagate(times, state0, 250.0, use_drag=False, backend="rust")
    assert np.allclose(py_states, rs_states, rtol=1e-11, atol=1e-11)
