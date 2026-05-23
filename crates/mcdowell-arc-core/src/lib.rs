use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

const MU_EARTH_KM3_S2: f64 = 398_600.4418;
const R_EARTH_KM: f64 = 6_378.137;
const OMEGA_EARTH_RAD_S: f64 = 7.292_115_9e-5;

const DENSITY_TABLE: &[(f64, f64)] = &[
    (0.0, 1.225),
    (25.0, 3.899e-2),
    (30.0, 1.774e-2),
    (40.0, 3.972e-3),
    (50.0, 1.057e-3),
    (60.0, 3.206e-4),
    (70.0, 8.770e-5),
    (80.0, 1.905e-5),
    (90.0, 3.396e-6),
    (100.0, 5.297e-7),
    (110.0, 9.661e-8),
    (120.0, 2.438e-8),
    (130.0, 8.484e-9),
    (140.0, 3.845e-9),
    (150.0, 2.070e-9),
    (180.0, 5.464e-10),
    (200.0, 2.789e-10),
    (250.0, 7.248e-11),
    (300.0, 2.418e-11),
    (350.0, 9.518e-12),
    (400.0, 3.725e-12),
    (450.0, 1.585e-12),
    (500.0, 6.967e-13),
    (600.0, 1.454e-13),
    (700.0, 3.614e-14),
    (800.0, 1.170e-14),
    (900.0, 5.245e-15),
    (1000.0, 3.019e-15),
];

#[inline]
fn check_finite(value: f64, name: &str) -> PyResult<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(PyValueError::new_err(format!("{name} must be finite")))
    }
}

#[inline]
fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec_add(a: &[f64; 6], b: &[f64; 6], scale: f64) -> [f64; 6] {
    [
        a[0] + scale * b[0],
        a[1] + scale * b[1],
        a[2] + scale * b[2],
        a[3] + scale * b[3],
        a[4] + scale * b[4],
        a[5] + scale * b[5],
    ]
}

fn density_scalar(altitude_km: f64) -> f64 {
    let mut alt = altitude_km;
    let min_alt = DENSITY_TABLE[0].0;
    let max_alt = DENSITY_TABLE[DENSITY_TABLE.len() - 1].0;
    if alt < min_alt {
        alt = min_alt;
    }
    if alt > max_alt {
        alt = max_alt;
    }

    for window in DENSITY_TABLE.windows(2) {
        let (h0, rho0) = window[0];
        let (h1, rho1) = window[1];
        if alt >= h0 && alt <= h1 {
            let frac = if h1 > h0 { (alt - h0) / (h1 - h0) } else { 0.0 };
            let log_rho = rho0.ln() + frac * (rho1.ln() - rho0.ln());
            return log_rho.exp();
        }
    }
    DENSITY_TABLE[DENSITY_TABLE.len() - 1].1
}

fn rotation_velocity(r: &[f64; 3]) -> [f64; 3] {
    [-OMEGA_EARTH_RAD_S * r[1], OMEGA_EARTH_RAD_S * r[0], 0.0]
}

fn acceleration(state: &[f64; 6], ballistic_coef_kg_m2: f64, use_drag: bool) -> PyResult<[f64; 3]> {
    let r = [state[0], state[1], state[2]];
    let v = [state[3], state[4], state[5]];
    let r_norm = norm3(&r);
    if r_norm <= 0.0 || !r_norm.is_finite() {
        return Err(PyValueError::new_err("position norm must be positive and finite"));
    }

    let mut a = [
        -MU_EARTH_KM3_S2 * r[0] / r_norm.powi(3),
        -MU_EARTH_KM3_S2 * r[1] / r_norm.powi(3),
        -MU_EARTH_KM3_S2 * r[2] / r_norm.powi(3),
    ];

    if !use_drag {
        return Ok(a);
    }
    if ballistic_coef_kg_m2 <= 0.0 || !ballistic_coef_kg_m2.is_finite() {
        return Ok(a);
    }

    let alt_km = (r_norm - R_EARTH_KM).max(0.0);
    let rho = density_scalar(alt_km);
    let v_atm = rotation_velocity(&r);
    let v_rel = [v[0] - v_atm[0], v[1] - v_atm[1], v[2] - v_atm[2]];
    let v_rel_norm_km_s = norm3(&v_rel);
    let v_rel_m_s = v_rel_norm_km_s * 1000.0;

    if v_rel_m_s <= 0.0 || !v_rel_m_s.is_finite() {
        return Ok(a);
    }

    let drag_mag_m_s2 = -0.5 * rho * v_rel_m_s * v_rel_m_s / ballistic_coef_kg_m2;
    let unit = [v_rel[0] / v_rel_norm_km_s, v_rel[1] / v_rel_norm_km_s, v_rel[2] / v_rel_norm_km_s];
    a[0] += drag_mag_m_s2 * unit[0] / 1000.0;
    a[1] += drag_mag_m_s2 * unit[1] / 1000.0;
    a[2] += drag_mag_m_s2 * unit[2] / 1000.0;
    Ok(a)
}

fn rhs(state: &[f64; 6], ballistic_coef_kg_m2: f64, use_drag: bool) -> PyResult<[f64; 6]> {
    let a = acceleration(state, ballistic_coef_kg_m2, use_drag)?;
    Ok([state[3], state[4], state[5], a[0], a[1], a[2]])
}

fn rk4_step(state: &[f64; 6], dt_s: f64, ballistic_coef_kg_m2: f64, use_drag: bool) -> PyResult<[f64; 6]> {
    let k1 = rhs(state, ballistic_coef_kg_m2, use_drag)?;
    let s2 = vec_add(state, &k1, 0.5 * dt_s);
    let k2 = rhs(&s2, ballistic_coef_kg_m2, use_drag)?;
    let s3 = vec_add(state, &k2, 0.5 * dt_s);
    let k3 = rhs(&s3, ballistic_coef_kg_m2, use_drag)?;
    let s4 = vec_add(state, &k3, dt_s);
    let k4 = rhs(&s4, ballistic_coef_kg_m2, use_drag)?;

    let mut out = [0.0; 6];
    for i in 0..6 {
        out[i] = state[i] + (dt_s / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    Ok(out)
}

fn state6_from_vec(values: Vec<f64>, name: &str) -> PyResult<[f64; 6]> {
    if values.len() != 6 {
        return Err(PyValueError::new_err(format!("{name} must contain exactly 6 values")));
    }
    for (idx, value) in values.iter().enumerate() {
        check_finite(*value, &format!("{name}[{idx}]"))?;
    }
    Ok([values[0], values[1], values[2], values[3], values[4], values[5]])
}

#[pyfunction]
fn density_kg_m3(altitude_km: f64) -> PyResult<f64> {
    check_finite(altitude_km, "altitude_km")?;
    Ok(density_scalar(altitude_km))
}

#[pyfunction]
fn earth_rotation_velocity_km_s(r_eci_km: Vec<f64>) -> PyResult<Vec<f64>> {
    if r_eci_km.len() != 3 {
        return Err(PyValueError::new_err("r_eci_km must contain exactly 3 values"));
    }
    for (idx, value) in r_eci_km.iter().enumerate() {
        check_finite(*value, &format!("r_eci_km[{idx}]"))?;
    }
    let r = [r_eci_km[0], r_eci_km[1], r_eci_km[2]];
    Ok(rotation_velocity(&r).to_vec())
}

#[pyfunction]
fn propagate(
    sample_times_s: Vec<f64>,
    state0: Vec<f64>,
    ballistic_coef_kg_m2: f64,
    use_drag: bool,
    max_step_s: f64,
) -> PyResult<Vec<Vec<f64>>> {
    if sample_times_s.is_empty() {
        return Err(PyValueError::new_err("sample_times_s must be non-empty"));
    }
    if max_step_s <= 0.0 || !max_step_s.is_finite() {
        return Err(PyValueError::new_err("max_step_s must be positive and finite"));
    }
    check_finite(ballistic_coef_kg_m2, "ballistic_coef_kg_m2")?;
    for (idx, t) in sample_times_s.iter().enumerate() {
        check_finite(*t, &format!("sample_times_s[{idx}]"))?;
        if idx > 0 && *t < sample_times_s[idx - 1] {
            return Err(PyValueError::new_err("sample_times_s must be sorted ascending"));
        }
    }

    let mut state = state6_from_vec(state0, "state0")?;
    let mut output = Vec::with_capacity(sample_times_s.len());
    output.push(state.to_vec());
    let mut t_current = sample_times_s[0];

    for t_target in sample_times_s.iter().skip(1) {
        while t_current < *t_target - 1e-12 {
            let dt = max_step_s.min(*t_target - t_current);
            state = rk4_step(&state, dt, ballistic_coef_kg_m2, use_drag)?;
            t_current += dt;
            if state.iter().any(|v| !v.is_finite()) {
                return Err(PyRuntimeError::new_err("propagation produced a non-finite state"));
            }
            let r_norm = norm3(&[state[0], state[1], state[2]]);
            if r_norm < R_EARTH_KM + 50.0 {
                return Err(PyRuntimeError::new_err("propagation hit the reentry guard before all sample times"));
            }
        }
        output.push(state.to_vec());
    }

    Ok(output)
}

#[pyfunction]
fn summarize_orbit(state: Vec<f64>) -> PyResult<(f64, Option<f64>, f64, Option<f64>, f64)> {
    let state = state6_from_vec(state, "state")?;
    let r = [state[0], state[1], state[2]];
    let v = [state[3], state[4], state[5]];
    let r_norm = norm3(&r);
    let v_norm = norm3(&v);
    let h_vec = cross(&r, &v);
    let h = norm3(&h_vec);
    if r_norm <= 0.0 || h <= 0.0 {
        return Err(PyValueError::new_err("degenerate orbit state"));
    }

    let energy = 0.5 * v_norm * v_norm - MU_EARTH_KM3_S2 / r_norm;
    let vxh = cross(&v, &h_vec);
    let e_vec = [
        vxh[0] / MU_EARTH_KM3_S2 - r[0] / r_norm,
        vxh[1] / MU_EARTH_KM3_S2 - r[1] / r_norm,
        vxh[2] / MU_EARTH_KM3_S2 - r[2] / r_norm,
    ];
    let ecc = norm3(&e_vec);
    let p = h * h / MU_EARTH_KM3_S2;
    let rp = p / (1.0 + ecc);
    let perigee_km = rp - R_EARTH_KM;

    let (apogee_km, semi_major_axis_km) = if energy < 0.0 {
        let a = -MU_EARTH_KM3_S2 / (2.0 * energy);
        let ra = a * (1.0 + ecc);
        (Some(ra - R_EARTH_KM), Some(a))
    } else {
        (None, None)
    };

    Ok((perigee_km, apogee_km, ecc, semi_major_axis_km, energy))
}

#[pymodule]
fn mcdowell_arc_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("MU_EARTH_KM3_S2", MU_EARTH_KM3_S2)?;
    m.add("R_EARTH_KM", R_EARTH_KM)?;
    m.add("OMEGA_EARTH_RAD_S", OMEGA_EARTH_RAD_S)?;
    m.add_function(wrap_pyfunction!(density_kg_m3, m)?)?;
    m.add_function(wrap_pyfunction!(earth_rotation_velocity_km_s, m)?)?;
    m.add_function(wrap_pyfunction!(propagate, m)?)?;
    m.add_function(wrap_pyfunction!(summarize_orbit, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn density_decreases_with_altitude() {
        assert!(density_scalar(100.0) > density_scalar(150.0));
        assert!(density_scalar(150.0) > density_scalar(250.0));
    }

    #[test]
    fn one_step_propagation_returns_two_rows() {
        let r = R_EARTH_KM + 400.0;
        let v = (MU_EARTH_KM3_S2 / r).sqrt();
        let out = propagate(vec![0.0, 10.0], vec![r, 0.0, 0.0, 0.0, v, 0.0], 250.0, false, 1.0)
            .expect("propagation should succeed");
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].len(), 6);
    }
}
