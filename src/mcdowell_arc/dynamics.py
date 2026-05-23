"""Orbital dynamics with a drag term in a co-rotating atmosphere."""

from __future__ import annotations

import numpy as np

from .atmosphere import density_kg_m3
from .constants import MU_EARTH_KM3_S2, R_EARTH_KM
from .fast import propagate_core, resolve_backend, validate_backend
from .frames import earth_rotation_velocity_km_s


def acceleration_km_s2(
    t_s: float,
    state: np.ndarray,
    ballistic_coef_kg_m2: float,
    use_drag: bool = True,
) -> np.ndarray:
    """Compute gravity plus drag acceleration for a 6-state vector.

    `ballistic_coef_kg_m2` is beta = mass / (Cd * area). Larger beta means less drag.
    """
    del t_s
    state = np.asarray(state, dtype=float)
    if state.shape != (6,):
        raise ValueError("state must be a 6-element vector [x,y,z,vx,vy,vz].")
    r = state[:3]
    v = state[3:]
    r_norm = np.linalg.norm(r)
    if r_norm <= 0 or not np.isfinite(r_norm):
        raise ValueError("Position norm must be positive and finite.")

    a_gravity = -MU_EARTH_KM3_S2 * r / r_norm**3
    if not use_drag:
        return a_gravity

    alt_km = max(r_norm - R_EARTH_KM, 0.0)
    rho = density_kg_m3(alt_km)
    v_atm = earth_rotation_velocity_km_s(r)
    v_rel_km_s = v - v_atm
    v_rel_norm_km_s = np.linalg.norm(v_rel_km_s)
    v_rel_m_s = v_rel_norm_km_s * 1000.0

    if v_rel_m_s <= 0.0 or ballistic_coef_kg_m2 <= 0.0:
        return a_gravity

    # Drag acceleration magnitude: 0.5 * rho * v^2 / beta in m/s^2.
    a_drag_m_s2_vec = -0.5 * rho * v_rel_m_s**2 / ballistic_coef_kg_m2 * (
        v_rel_km_s / v_rel_norm_km_s
    )
    a_drag_km_s2 = a_drag_m_s2_vec / 1000.0
    return a_gravity + a_drag_km_s2


def rhs(t_s: float, state: np.ndarray, ballistic_coef_kg_m2: float, use_drag: bool = True) -> np.ndarray:
    """First-order ODE right-hand side."""
    return np.hstack([state[3:], acceleration_km_s2(t_s, state, ballistic_coef_kg_m2, use_drag=use_drag)])


def _rk4_step(t_s: float, state: np.ndarray, dt_s: float, ballistic_coef_kg_m2: float, use_drag: bool) -> np.ndarray:
    """One fixed RK4 step."""
    k1 = rhs(t_s, state, ballistic_coef_kg_m2, use_drag=use_drag)
    k2 = rhs(t_s + 0.5 * dt_s, state + 0.5 * dt_s * k1, ballistic_coef_kg_m2, use_drag=use_drag)
    k3 = rhs(t_s + 0.5 * dt_s, state + 0.5 * dt_s * k2, ballistic_coef_kg_m2, use_drag=use_drag)
    k4 = rhs(t_s + dt_s, state + dt_s * k3, ballistic_coef_kg_m2, use_drag=use_drag)
    return state + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _propagate_python(
    sample_times_s: np.ndarray,
    state0: np.ndarray,
    ballistic_coef_kg_m2: float,
    use_drag: bool = True,
    max_step_s: float = 1.0,
) -> np.ndarray:
    """Pure-Python/NumPy propagation implementation."""
    times = np.asarray(sample_times_s, dtype=float)
    if times.ndim != 1 or len(times) == 0:
        raise ValueError("sample_times_s must be a non-empty 1-D array.")
    if not np.all(np.isfinite(times)):
        raise ValueError("sample_times_s must contain finite values.")
    if np.any(np.diff(times) < 0):
        raise ValueError("sample_times_s must be sorted ascending.")
    if max_step_s <= 0 or not np.isfinite(max_step_s):
        raise ValueError("max_step_s must be positive and finite.")

    state = np.asarray(state0, dtype=float).copy()
    if state.shape != (6,):
        raise ValueError("state0 must be a 6-element vector [x,y,z,vx,vy,vz].")
    if not np.all(np.isfinite(state)):
        raise ValueError("state0 must contain only finite values.")

    out = [state.copy()]
    t_current = float(times[0])

    for t_target in times[1:]:
        t_target = float(t_target)
        while t_current < t_target - 1e-12:
            dt = min(max_step_s, t_target - t_current)
            state = _rk4_step(t_current, state, dt, ballistic_coef_kg_m2, use_drag)
            t_current += dt
            if not np.all(np.isfinite(state)):
                raise RuntimeError("Propagation produced a non-finite state.")
            if np.linalg.norm(state[:3]) < R_EARTH_KM + 50.0:
                raise RuntimeError("Propagation hit the reentry guard before all sample times.")
        out.append(state.copy())

    return np.vstack(out)


def propagate(
    sample_times_s: np.ndarray,
    state0: np.ndarray,
    ballistic_coef_kg_m2: float,
    use_drag: bool = True,
    max_step_s: float = 1.0,
    backend: str = "auto",
) -> np.ndarray:
    """Propagate state0 to the requested sample times using fixed-step RK4.

    `backend="auto"` uses the Rust core when installed and falls back to the
    Python implementation otherwise. `backend="rust"` fails clearly when the
    extension is unavailable.
    """
    selected = validate_backend(backend)
    if resolve_backend(selected) == "rust":
        return propagate_core(sample_times_s, state0, ballistic_coef_kg_m2, use_drag, max_step_s)
    return _propagate_python(sample_times_s, state0, ballistic_coef_kg_m2, use_drag, max_step_s)
