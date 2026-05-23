# Scientific Model Notes

## Current model

The current release uses a reduced two-dimensional ECI state:

```text
state0 = [x, y, z, vx, vy, vz]
       = [R_earth + altitude0 + offset, 0, 0, vr, vt, 0]
```

The observations are:

- elapsed time in seconds
- altitude in kilometers
- Earth-relative speed magnitude in kilometers per second

The model integrates gravity plus optional drag in a co-rotating atmosphere.

## Why this is still not publication-grade

Altitude and speed magnitude alone do not uniquely determine a full 3-D inertial state. A serious trajectory fit needs additional geometry such as:

- latitude and longitude
- range and range rate
- azimuth/elevation from a known observer
- launch site and launch azimuth assumptions
- inclination or plane constraints
- external tracking observations

The current model is useful for software plumbing, sensitivity testing, public telemetry experiments, and uncertainty framing. It should not be represented as a definitive orbital reconstruction.

## Atmosphere caveat

The built-in atmosphere table is a rough log-interpolated density table. It should be replaced by a vetted model before high-confidence scientific claims.

## Correct public framing

Strong framing:

> This is an uncertainty-aware experimental fit under a reduced observation model.

Weak or misleading framing:

> This proves the object was orbital/non-orbital.
