# Validation Plan

This project should earn scientific credibility through layered validation.

## Level 0: software sanity

- unit tests for frame helpers
- unit tests for orbit summaries
- Python/Rust propagation agreement tests
- CLI smoke tests
- CSV validation tests

## Level 1: synthetic truth

Generate synthetic states with known perigee/apogee, then sample altitude/speed with noise. The fitter should recover the known state within uncertainty bounds under the same model.

## Level 2: model sensitivity

Run sensitivity sweeps over:

- ballistic coefficient
- drag on/off
- time-window selection
- altitude uncertainty
- speed uncertainty
- atmosphere threshold

## Level 3: independent event checks

Compare against public cases where an independent orbit estimate exists. Report assumptions and residuals, not just final perigee.

## Level 4: full geometry

Move beyond webcast altitude/speed by supporting real observation geometry and Earth orientation handling.
