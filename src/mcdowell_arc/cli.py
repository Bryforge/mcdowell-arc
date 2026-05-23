"""Command-line interface for McDowell Arc."""

from __future__ import annotations

import argparse
import json
import platform
import sys

import numpy as np
import pandas as pd
import scipy

from . import __version__
from .fast import core_available, core_version, resolve_backend
from .fit import fit_observations, load_observations_csv, monte_carlo


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mcdowell-arc")
    parser.add_argument("--version", action="version", version=f"mcdowell-arc {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    fit = sub.add_parser("fit", help="fit a drag-aware ECI trajectory from webcast CSV")
    fit.add_argument("csv", help="CSV containing t_s, altitude_km, speed_km_s")
    fit.add_argument("--monte-carlo", type=int, default=0, help="number of Monte Carlo refits")
    fit.add_argument("--atmosphere-km", type=float, default=100.0, help="perigee threshold for atmosphere crossing")
    fit.add_argument("--no-drag", action="store_true", help="disable drag for comparison")
    fit.add_argument("--seed", type=int, default=7, help="Monte Carlo random seed")
    fit.add_argument("--backend", choices=["auto", "python", "rust"], default="auto", help="propagation backend")
    fit.add_argument("--max-step-s", type=float, default=1.0, help="maximum fixed RK4 propagation step in seconds")
    fit.add_argument("--max-nfev", type=int, default=35, help="maximum optimizer function evaluations")

    sub.add_parser("doctor", help="show environment and backend availability")
    return parser


def _doctor_payload() -> dict:
    return {
        "mcdowell_arc_version": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "rust_core_available": core_available(),
        "rust_core_version": core_version(),
        "default_backend": resolve_backend("auto"),
    }


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "doctor":
        print(json.dumps(_doctor_payload(), indent=2, sort_keys=True))
        return 0

    if args.command == "fit":
        obs = load_observations_csv(args.csv)
        backend = resolve_backend(args.backend)
        nominal = fit_observations(
            obs,
            use_drag=not args.no_drag,
            backend=args.backend,
            max_step_s=args.max_step_s,
            max_nfev=args.max_nfev,
        )
        output = {
            "model": "reduced_2d_eci_drag" if not args.no_drag else "reduced_2d_eci_no_drag",
            "model_status": "engineering_mvp_not_publication_grade",
            "backend": backend,
            "nominal": {
                "perigee_km": nominal.orbit.perigee_km,
                "apogee_km": nominal.orbit.apogee_km,
                "eccentricity": nominal.orbit.eccentricity,
                "semi_major_axis_km": nominal.orbit.semi_major_axis_km,
                "ballistic_coef_kg_m2": nominal.ballistic_coef_kg_m2,
                "residual_rms_sigma": nominal.residual_rms,
                "optimizer_cost": nominal.optimizer_cost,
                "optimizer_success": nominal.optimizer_success,
                "optimizer_message": nominal.message,
            },
            "state0_eci_km_km_s": nominal.state0.tolist(),
            "interpretation": (
                "Report probability_perigee_above_threshold, not only the nominal perigee. "
                "This model remains reduced 2-D and needs full tracking geometry, atmosphere validation, "
                "and independent checks before publication-grade trajectory claims."
            ),
        }
        if args.monte_carlo > 0:
            output["monte_carlo"] = monte_carlo(
                obs,
                samples=args.monte_carlo,
                atmosphere_threshold_km=args.atmosphere_km,
                seed=args.seed,
                use_drag=not args.no_drag,
                backend=args.backend,
                max_step_s=args.max_step_s,
            )
        print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
