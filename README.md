# McDowell Arc

**McDowell Arc** is a small starter tool for drag-aware, Earth-centered inertial trajectory fitting from webcast-style telemetry. It is named as a working-title tribute to Jonathan McDowell's challenge: do better than a clean ballistic ellipse fitted near apogee.

The MVP takes timestamped altitude and Earth-relative speed samples, fits a reduced 2-D ECI trajectory with atmospheric drag, and reports a nominal perigee/apogee plus Monte Carlo confidence bounds. The primary result is not just “the fitted perigee,” but the probability that perigee rises above a configurable atmosphere threshold.

> Public naming note: before publishing under a living person's name, consider asking permission or using a neutral package name with an acknowledgement section.

## What it answers

Instead of asking only:

> What ellipse fits these samples?

McDowell Arc asks:

> Given measurement uncertainty and drag, what is the probability the fitted perigee remains above the atmosphere?

## Current MVP scope

Implemented now:

- Webcast CSV ingestion: `t_s`, `altitude_km`, `speed_km_s`
- Optional uncertainty columns: `sigma_altitude_km`, `sigma_speed_km_s`
- Earth rotation correction in the drag-relative velocity term
- Simplified co-rotating atmosphere
- Built-in rough density table for 0-1000 km altitude
- Weighted least-squares fitting
- Monte Carlo uncertainty propagation
- Perigee/apogee reporting

Important limitations:

- The current fitter is a reduced 2-D approximation. Full ECI reconstruction needs position/velocity direction information, such as lat/lon/range/range-rate, tracking geometry, or assumptions about launch azimuth/inclination.
- The built-in atmosphere is intentionally simple. Replace it with NRLMSISE-00 or another vetted density model before serious public claims.
- Near-apogee-only telemetry can be underconstrained. The tool reports uncertainty so that weak fits remain visibly weak.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Run the example

```bash
mcdowell-arc fit examples/sample_webcast.csv --monte-carlo 20 --atmosphere-km 100
```

## CSV format

Required columns:

```csv
t_s,altitude_km,speed_km_s
0,210,7.81
10,211,7.80
20,210.8,7.79
```

Recommended uncertainty columns:

```csv
t_s,altitude_km,speed_km_s,sigma_altitude_km,sigma_speed_km_s
0,210,7.81,3,0.05
10,211,7.80,3,0.05
```

## Roadmap

1. Add full ECEF/geodetic to ECI observation support.
2. Add real Earth orientation / sidereal-time handling.
3. Add NRLMSISE-00 density as a plug-in backend.
4. Support boost/coast/drag phase boundaries.
5. Add covariance estimation from the least-squares Jacobian.
6. Add plots for altitude, speed, residuals, and perigee distribution.
7. Add a reproducible report generator suitable for public trajectory notes.

