"""Optional Rust acceleration bridge for McDowell Arc.

The package is fully usable without Rust. When the `mcdowell_arc_core` extension
is installed, selected deterministic kernels can be delegated to it.
"""

from __future__ import annotations

from functools import lru_cache
from importlib import import_module
from types import ModuleType
from typing import Literal

import numpy as np

BackendName = Literal["auto", "python", "rust"]


class BackendUnavailableError(RuntimeError):
    """Raised when the Rust backend is explicitly requested but unavailable."""


@lru_cache(maxsize=1)
def _load_core() -> ModuleType | None:
    try:
        return import_module("mcdowell_arc_core")
    except ModuleNotFoundError:
        return None


def core_available() -> bool:
    """Return True when the optional Rust extension can be imported."""
    return _load_core() is not None


def core_version() -> str | None:
    """Return the Rust core version, if the extension is installed."""
    core = _load_core()
    return None if core is None else str(getattr(core, "__version__", "unknown"))


def require_core() -> ModuleType:
    """Return the Rust extension or raise a clear installation error."""
    core = _load_core()
    if core is None:
        raise BackendUnavailableError(
            "Rust backend requested, but mcdowell_arc_core is not installed. "
            "Run: maturin develop --manifest-path crates/mcdowell-arc-core/Cargo.toml --release"
        )
    return core


def validate_backend(backend: str) -> BackendName:
    """Validate a backend selector from CLI/API input."""
    allowed = {"auto", "python", "rust"}
    if backend not in allowed:
        raise ValueError(f"backend must be one of {sorted(allowed)}")
    return backend  # type: ignore[return-value]


def resolve_backend(backend: str) -> Literal["python", "rust"]:
    """Resolve `auto` to the concrete backend that will be used."""
    selected = validate_backend(backend)
    if selected == "python":
        return "python"
    if selected == "rust":
        require_core()
        return "rust"
    return "rust" if core_available() else "python"


def propagate_core(
    sample_times_s: np.ndarray,
    state0: np.ndarray,
    ballistic_coef_kg_m2: float,
    use_drag: bool,
    max_step_s: float,
) -> np.ndarray:
    """Call Rust propagation and return a NumPy array."""
    core = require_core()
    rows = core.propagate(
        np.asarray(sample_times_s, dtype=float).tolist(),
        np.asarray(state0, dtype=float).tolist(),
        float(ballistic_coef_kg_m2),
        bool(use_drag),
        float(max_step_s),
    )
    return np.asarray(rows, dtype=float)
