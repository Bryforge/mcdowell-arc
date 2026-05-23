"""Command-line interface for McDowell Arc."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from .fit import fit_observations, load_observations_csv, monte_carlo


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mcdowell-arc")
    sub = parser.add_subparsers(dest="command", required=True)

    fit = sub.add_parser("fit", help="fit a drag-aware ECI trajectory from webcast CSV")
    fit.add_argument("csv", help="CSV containing t_s, altitude_km, speed_km_s")
    fit.add_argument("--monte-carlo", type=int, default=0, help="number of Monte Carlo refits")
    fit.add_argument("--atmosphere-km", type=float, default=100.0, help="perigee threshold for atmosphere crossing")
    fit.add_argument("--no-drag", action="store_true", help="disable drag for comparison")
    fit.add_argument("--seed", type=int, default=7, help="Monte Carlo random seed")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "fit":
        obs = load_observations_csv(args.csv)
        nominal = fit_observations(obs, use_drag=not args.no_drag)
        output = {
            "model": "reduced_2d_eci_drag_mvp" if not args.no_drag else "reduced_2d_eci_no_drag_mvp",
            "nominal": {
                "perigee_km": nominal.orbit.perigee_km,
                "apogee_km": nominal.orbit.apogee_km,
                "eccentricity": nominal.orbit.eccentricity,
                "semi_major_axis_km": nominal.orbit.semi_major_axis_km,
                "ballistic_coef_kg_m2": nominal.ballistic_coef_kg_m2,
                "residual_rms_sigma": nominal.residual_rms,
                "optimizer_success": nominal.optimizer_success,
                "optimizer_message": nominal.message,
            },
            "state0_eci_km_km_s": nominal.state0.tolist(),
            "interpretation": (
                "Report probability_perigee_above_threshold, not only the nominal perigee. "
                "MVP is reduced 2-D and needs full tracking geometry for publication-grade results."
            ),
        }
        if args.monte_carlo > 0:
            output["monte_carlo"] = monte_carlo(
                obs,
                samples=args.monte_carlo,
                atmosphere_threshold_km=args.atmosphere_km,
                seed=args.seed,
                use_drag=not args.no_drag,
            )
        print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
