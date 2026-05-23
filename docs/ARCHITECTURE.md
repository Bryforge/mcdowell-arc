# Architecture

McDowell Arc is split into two layers.

## Python layer: `mcdowell_arc`

Responsibilities:

- command-line interface
- CSV ingestion and validation
- weighted least-squares fitting
- Monte Carlo uncertainty orchestration
- JSON reporting
- plotting and report generation in future releases

## Rust layer: `mcdowell-arc-core`

Responsibilities:

- deterministic propagation kernels
- drag acceleration math
- orbit summary math
- future parallel Monte Carlo kernels
- future optimized observation transforms

The Python package remains usable without Rust. The CLI option `--backend auto` selects Rust when the extension is installed and Python otherwise.

## Why this split exists

Python is better for fast scientific iteration and reproducible user workflows. Rust is better for deterministic, tested numerical kernels that may need to run thousands or millions of times. The design goal is not to replace Python, but to move hot loops into a small, auditable core.
