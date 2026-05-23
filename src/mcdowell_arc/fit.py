"""Weighted fitting and uncertainty propagation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .constants import R_EARTH_KM
from .dynamics import propagate
from .fast import validate_backend
from .frames import earth_rotation_velocity_km_s
from .orbit import OrbitSummary, summarize_orbit


@dataclass(frozen=True)
class Observations:
    """Webcast-style trajectory observations.

    The MVP observation format intentionally contains only time, altitude, and
    Earth-relative speed. That makes it useful for public webcast telemetry, but
    it is underconstrained compared with real tracking geometry.
    """

    t_s: np.ndarray
    altitude_km: np.ndarray
    speed_km_s: np.ndarray
    sigma_altitude_km: np.ndarray
    sigma_speed_km_s: np.ndarray

    def validated(self) -> Observations:
        arrays = {
            "t_s": np.asarray(self.t_s, dtype=float),
            "altitude_km": np.asarray(self.altitude_km, dtype=float),
            "speed_km_s": np.asarray(self.speed_km_s, dtype=float),
            "sigma_altitude_km": np.asarray(self.sigma_altitude_km, dtype=float),
            "sigma_speed_km_s": np.asarray(self.sigma_speed_km_s, dtype=float),
        }
        length = len(arrays["t_s"])
        if length < 3:
            raise ValueError("At least three observations are required for fitting.")
        for name, values in arrays.items():
            if values.ndim != 1:
                raise ValueError(f"{name} must be a 1-D array.")
            if len(values) != length:
                raise ValueError(f"{name} length must match t_s length.")
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must contain only finite values.")
        if np.any(np.diff(arrays["t_s"]) <= 0):
            raise ValueError("t_s must be strictly increasing after sorting.")
        if np.any(arrays["altitude_km"] < -1.0):
            raise ValueError("altitude_km contains values below the accepted guard band.")
        if np.any(arrays["speed_km_s"] <= 0.0):
            raise ValueError("speed_km_s must be positive.")
        if np.any(arrays["sigma_altitude_km"] <= 0.0) or np.any(arrays["sigma_speed_km_s"] <= 0.0):
            raise ValueError("uncertainty columns must be positive.")
        return Observations(**arrays)


@dataclass(frozen=True)
class FitOptions:
    """Controls for deterministic fitting and propagation."""

    use_drag: bool = True
    backend: str = "auto"
    max_step_s: float = 1.0
    max_nfev: int = 35

    def validated(self) -> FitOptions:
        validate_backend(self.backend)
        if self.max_step_s <= 0 or not np.isfinite(self.max_step_s):
            raise ValueError("max_step_s must be positive and finite.")
        if self.max_nfev <= 0:
            raise ValueError("max_nfev must be positive.")
        return self


@dataclass(frozen=True)
class FitResult:
    state0: np.ndarray
    ballistic_coef_kg_m2: float
    orbit: OrbitSummary
    residual_rms: float
    optimizer_cost: float
    optimizer_success: bool
    message: str


def _numeric_series(df: pd.DataFrame, column: str) -> np.ndarray:
    try:
        return pd.to_numeric(df[column], errors="raise").to_numpy(float)
    except Exception as exc:
        raise ValueError(f"CSV column {column!r} must be numeric.") from exc


def load_observations_csv(path: str | Path) -> Observations:
    """Load webcast-style observations from CSV and validate them."""
    df = pd.read_csv(path)
    required = {"t_s", "altitude_km", "speed_km_s"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required CSV columns: {', '.join(missing)}")

    df = df.sort_values("t_s").reset_index(drop=True)
    sigma_alt = _numeric_series(df, "sigma_altitude_km") if "sigma_altitude_km" in df else np.full(len(df), 3.0)
    sigma_speed = _numeric_series(df, "sigma_speed_km_s") if "sigma_speed_km_s" in df else np.full(len(df), 0.05)
    sigma_alt = np.where(sigma_alt > 0, sigma_alt, 3.0)
    sigma_speed = np.where(sigma_speed > 0, sigma_speed, 0.05)

    t = _numeric_series(df, "t_s")
    t = t - t[0]
    obs = Observations(
        t_s=t,
        altitude_km=_numeric_series(df, "altitude_km"),
        speed_km_s=_numeric_series(df, "speed_km_s"),
        sigma_altitude_km=sigma_alt,
        sigma_speed_km_s=sigma_speed,
    )
    return obs.validated()


def _params_to_state(params: np.ndarray, obs: Observations) -> tuple[np.ndarray, float]:
    """Map optimizer parameters into a 2-D initial ECI state.

    params = [r0_offset_km, vr0_km_s, vt0_km_s, log_beta]
    """
    r0_km = R_EARTH_KM + obs.altitude_km[0] + params[0]
    state0 = np.array([r0_km, 0.0, 0.0, params[1], params[2], 0.0], dtype=float)
    beta = float(np.exp(params[3]))
    return state0, beta


def observables_from_states(states: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute altitude and Earth-relative speed from propagated ECI states."""
    states = np.asarray(states, dtype=float)
    r = states[:, :3]
    v = states[:, 3:]
    altitude = np.linalg.norm(r, axis=1) - R_EARTH_KM
    v_atm = np.array([earth_rotation_velocity_km_s(ri) for ri in r])
    rel_speed = np.linalg.norm(v - v_atm, axis=1)
    return altitude, rel_speed


def residuals(params: np.ndarray, obs: Observations, options: FitOptions) -> np.ndarray:
    """Weighted residual vector for least-squares fitting."""
    state0, beta = _params_to_state(params, obs)
    try:
        states = propagate(
            obs.t_s,
            state0,
            beta,
            use_drag=options.use_drag,
            max_step_s=options.max_step_s,
            backend=options.backend,
        )
        alt_model, speed_model = observables_from_states(states)
    except Exception:
        # Penalize impossible trial states without crashing the optimizer.
        return np.full(obs.t_s.size * 2, 1e6)

    return np.concatenate(
        [
            (alt_model - obs.altitude_km) / obs.sigma_altitude_km,
            (speed_model - obs.speed_km_s) / obs.sigma_speed_km_s,
        ]
    )


def initial_guess(obs: Observations) -> np.ndarray:
    """Build a simple first guess near the first observation."""
    obs = obs.validated()
    if len(obs.t_s) >= 2 and obs.t_s[-1] > obs.t_s[0]:
        vr_guess = np.gradient(obs.altitude_km, obs.t_s)[0]
    else:
        vr_guess = 0.0

    # At the equator, the atmosphere co-rotates at roughly 0.465 km/s.
    # Add it back because the state is inertial, while observations are Earth-relative speeds.
    vt_guess = float(max(obs.speed_km_s[0] + 0.465, 0.01))
    beta_guess = 250.0
    return np.array([0.0, vr_guess, vt_guess, np.log(beta_guess)], dtype=float)


def fit_observations(
    obs: Observations,
    use_drag: bool = True,
    initial_params: np.ndarray | None = None,
    max_nfev: int = 35,
    backend: str = "auto",
    max_step_s: float = 1.0,
) -> FitResult:
    """Fit a reduced 2-D drag-aware ECI trajectory."""
    obs = obs.validated()
    options = FitOptions(use_drag=use_drag, backend=backend, max_step_s=max_step_s, max_nfev=max_nfev).validated()
    x0 = initial_guess(obs) if initial_params is None else np.asarray(initial_params, dtype=float)
    bounds = (
        np.array([-100.0, -5.0, 0.1, np.log(10.0)]),
        np.array([100.0, 5.0, 12.0, np.log(20000.0)]),
    )
    opt = least_squares(
        residuals,
        x0=x0,
        bounds=bounds,
        args=(obs, options),
        loss="soft_l1",
        f_scale=2.0,
        max_nfev=options.max_nfev,
        x_scale=np.array([10.0, 1.0, 1.0, 2.0]),
    )
    state0, beta = _params_to_state(opt.x, obs)
    orbit = summarize_orbit(state0)
    res = residuals(opt.x, obs, options)
    return FitResult(
        state0=state0,
        ballistic_coef_kg_m2=beta,
        orbit=orbit,
        residual_rms=float(np.sqrt(np.mean(res**2))),
        optimizer_cost=float(opt.cost),
        optimizer_success=bool(opt.success),
        message=str(opt.message),
    )


def monte_carlo(
    obs: Observations,
    samples: int,
    atmosphere_threshold_km: float = 100.0,
    seed: int = 7,
    use_drag: bool = True,
    backend: str = "auto",
    max_step_s: float = 1.0,
) -> dict:
    """Perturb observations by their sigmas and refit many times."""
    obs = obs.validated()
    if samples <= 0:
        raise ValueError("samples must be positive.")
    rng = np.random.default_rng(seed)
    nominal = fit_observations(obs, use_drag=use_drag, backend=backend, max_step_s=max_step_s)
    perigees: list[float] = []
    apogees: list[float] = []
    betas: list[float] = []
    failures = 0

    for _ in range(samples):
        perturbed = Observations(
            t_s=obs.t_s.copy(),
            altitude_km=obs.altitude_km + rng.normal(0.0, obs.sigma_altitude_km),
            speed_km_s=obs.speed_km_s + rng.normal(0.0, obs.sigma_speed_km_s),
            sigma_altitude_km=obs.sigma_altitude_km.copy(),
            sigma_speed_km_s=obs.sigma_speed_km_s.copy(),
        ).validated()
        try:
            warm_start = np.array(
                [
                    nominal.state0[0] - (R_EARTH_KM + perturbed.altitude_km[0]),
                    nominal.state0[3],
                    nominal.state0[4],
                    np.log(nominal.ballistic_coef_kg_m2),
                ]
            )
            fit = fit_observations(
                perturbed,
                use_drag=use_drag,
                initial_params=warm_start,
                max_nfev=18,
                backend=backend,
                max_step_s=max_step_s,
            )
            perigees.append(fit.orbit.perigee_km)
            if fit.orbit.apogee_km is not None:
                apogees.append(fit.orbit.apogee_km)
            betas.append(fit.ballistic_coef_kg_m2)
        except Exception:
            failures += 1

    perigee_arr = np.asarray(perigees, dtype=float)
    apogee_arr = np.asarray(apogees, dtype=float)
    beta_arr = np.asarray(betas, dtype=float)
    if len(perigee_arr) == 0:
        raise RuntimeError("All Monte Carlo fits failed.")

    def pct(arr: np.ndarray, q: float) -> float | None:
        if len(arr) == 0:
            return None
        return float(np.percentile(arr, q))

    return {
        "samples_requested": int(samples),
        "samples_succeeded": int(len(perigee_arr)),
        "samples_failed": int(failures),
        "probability_perigee_above_threshold": float(np.mean(perigee_arr >= atmosphere_threshold_km)),
        "atmosphere_threshold_km": float(atmosphere_threshold_km),
        "perigee_km_p05": pct(perigee_arr, 5),
        "perigee_km_p50": pct(perigee_arr, 50),
        "perigee_km_p95": pct(perigee_arr, 95),
        "apogee_km_p05": pct(apogee_arr, 5),
        "apogee_km_p50": pct(apogee_arr, 50),
        "apogee_km_p95": pct(apogee_arr, 95),
        "ballistic_coef_kg_m2_p50": pct(beta_arr, 50),
    }
