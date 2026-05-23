# McDowell Arc

**McDowell Arc** is an experimental scientific software tool for drag-aware trajectory fitting and uncertainty analysis from simplified webcast-style trajectory observations.

It estimates nominal perigee/apogee values, confidence intervals, and the probability that a fitted perigee remains above a selected atmospheric threshold.

The v0.2 architecture separates the project into two parts:

- **`mcdowell-arc`** — Python CLI, CSV loading, reporting, validation, and scientific workflow.
- **`mcdowell-arc-core`** — Rust acceleration core exposed to Python through a native extension.

The goal is to keep the tool easy to use from Python while moving performance-sensitive trajectory calculations into a faster compiled backend.

> Named in honor of Jonathan McDowell's challenge to improve near-orbital trajectory fitting. This project is not affiliated with or endorsed by Jonathan McDowell unless explicitly stated.

---

## Current Status

McDowell Arc is currently an engineering MVP moving toward production-quality scientific software.

It includes:

- Webcast-style CSV ingestion
- Weighted trajectory fitting
- Simplified drag-aware dynamics
- Earth-rotation-aware drag-relative velocity
- Monte Carlo uncertainty propagation
- Nominal perigee/apogee reporting
- Confidence interval reporting
- Probability estimate for perigee above an atmospheric threshold
- Python backend
- Rust acceleration backend through `mcdowell-arc-core`
- Backend selection with `--backend auto|python|rust`
- Environment diagnostics with `mcdowell-arc doctor`
- Test coverage for Python and Rust-backed behavior

Scientific caveat: the current model is still simplified. It is useful for development, experimentation, validation workflows, and early trajectory-analysis research, but it is not yet publication-grade orbital determination software.

---

## Repository Layout

```text
mcdowell-arc/
├── mcdowell_arc/                  # Python package and CLI
├── crates/
│   └── mcdowell-arc-core/         # Rust/PyO3 acceleration core
├── examples/                      # Example CSV input data
├── tests/                         # Python test suite
├── docs/                          # Scientific and validation notes
├── pyproject.toml                 # Python package metadata
├── Cargo.toml                     # Rust workspace metadata
└── README.md
```

---

## Requirements

### Python

Python 3.10 or newer is recommended.

The project has been tested on Fedora with Python 3.14.

### Rust

The Rust backend requires Rust and Cargo.

On Fedora:

```bash
sudo dnf install -y rust cargo
```

You can also install Rust through `rustup`.

---

## Quick Start

Clone the repository:

```bash
git clone https://github.com/Bryforge/mcdowell-arc.git
cd mcdowell-arc
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Upgrade packaging tools:

```bash
python -m pip install --upgrade pip setuptools wheel
```

Install the Python package in editable development mode:

```bash
python -m pip install -e ".[dev]"
```

Run the tests:

```bash
python -m pytest -q
```

Run diagnostics:

```bash
mcdowell-arc doctor
```

At this point, the Python backend should work.

---

## Building the Rust Acceleration Core

Install `maturin`:

```bash
python -m pip install maturin
```

Build and install the Rust extension from inside the Rust crate directory:

```bash
cd crates/mcdowell-arc-core
maturin develop --release
cd ../..
```

Verify that both packages are installed:

```bash
python -m pip list | grep mcdowell
```

Expected result:

```text
mcdowell-arc       0.2.0
mcdowell-arc-core  0.2.0
```

Verify that Python can import both modules:

```bash
python -c "import mcdowell_arc; print('Python package loaded')"
python -c "import mcdowell_arc_core; print('Rust core loaded')"
```

Run the full test suite:

```bash
python -m pytest -q
```

Expected result:

```text
12 passed
```

---

## Usage

Run a fit using automatic backend selection:

```bash
mcdowell-arc fit examples/sample_webcast.csv --monte-carlo 1000 --atmosphere-km 100 --backend auto
```

Force the Python backend:

```bash
mcdowell-arc fit examples/sample_webcast.csv --monte-carlo 1000 --atmosphere-km 100 --backend python
```

Force the Rust backend:

```bash
mcdowell-arc fit examples/sample_webcast.csv --monte-carlo 1000 --atmosphere-km 100 --backend rust
```

If the Rust backend is not installed and `--backend rust` is selected, the command will fail. If `--backend auto` is selected, the tool will use Rust when available and fall back to Python when needed.

---

## Input CSV Format

The basic input CSV should contain:

```csv
t_s,altitude_km,speed_km_s
0,160.0,7.72
10,158.4,7.70
20,156.8,7.68
```

Required columns:

- `t_s` — time in seconds
- `altitude_km` — observed altitude in kilometers
- `speed_km_s` — observed speed in kilometers per second

Optional uncertainty columns may be added in future model expansions.

---

## Example Output

The CLI emits JSON output containing fields such as:

```json
{
  "optimizer_success": true,
  "backend": "rust",
  "samples_requested": 1000,
  "samples_succeeded": 1000,
  "samples_failed": 0,
  "nominal": {
    "perigee_km": 129.0,
    "apogee_km": 2242.0
  },
  "probability_perigee_above_threshold": 1.0
}
```

Exact values depend on the input data, random seed behavior, model settings, and backend.

---

## Development Workflow

Install development dependencies:

```bash
python -m pip install -e ".[dev]"
python -m pip install maturin
```

Build the Rust backend:

```bash
cd crates/mcdowell-arc-core
maturin develop --release
cd ../..
```

Run tests:

```bash
python -m pytest -q
```

Run diagnostics:

```bash
mcdowell-arc doctor
```

Run the example fit:

```bash
mcdowell-arc fit examples/sample_webcast.csv --monte-carlo 1000 --atmosphere-km 100 --backend auto
```

---

## Scientific Roadmap

Major next steps:

1. Full 3-D trajectory state modeling.
2. Observation geometry support: latitude, longitude, azimuth, range, range-rate, and launch-site assumptions.
3. Improved atmosphere modeling.
4. Better uncertainty modeling for webcast-derived measurements.
5. Validation against synthetic truth cases.
6. Validation against known historical trajectories where public data is available.
7. Parallel Monte Carlo execution.
8. Benchmarking Python backend versus Rust backend.
9. Plotting and report generation.
10. Formal model documentation.

---

## Production Readiness Goals

McDowell Arc should evolve toward:

- Reproducible scientific runs
- Versioned model assumptions
- Clear input schemas
- Deterministic test cases
- Synthetic validation datasets
- Backend parity tests
- CI testing for Python and Rust
- Structured JSON output
- Scientific documentation
- Honest uncertainty reporting

---

## License

Add a license before wider public release.
