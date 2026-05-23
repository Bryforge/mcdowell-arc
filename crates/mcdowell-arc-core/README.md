# mcdowell-arc-core

`mcdowell-arc-core` is the optional Rust acceleration layer for McDowell Arc.
It exposes a Python extension module named `mcdowell_arc_core` through PyO3.

The Python package works without this crate. Install the Rust core when you want
faster deterministic propagation inside the existing `mcdowell-arc` CLI.

```bash
python -m pip install maturin
maturin develop --manifest-path crates/mcdowell-arc-core/Cargo.toml --release
mcdowell-arc doctor
```
