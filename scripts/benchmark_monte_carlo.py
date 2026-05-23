from __future__ import annotations

import argparse
import time

from mcdowell_arc.fit import load_observations_csv, monte_carlo


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--backend", choices=["auto", "python", "rust"], default="auto")
    args = parser.parse_args()

    obs = load_observations_csv(args.csv)
    start = time.perf_counter()
    result = monte_carlo(obs, args.samples, backend=args.backend)
    elapsed = time.perf_counter() - start
    print(f"backend={args.backend} samples={args.samples} elapsed_s={elapsed:.3f}")
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
