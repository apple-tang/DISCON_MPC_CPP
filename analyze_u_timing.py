"""Summarize and plot MPC timing CSV files produced by the timing DLL."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


PHASE_COLUMNS = (
    "preprocess_us",
    "model_workspace_us",
    "prediction_us",
    "cost_us",
    "constraints_us",
    "qpoases_us",
    "extract_u_us",
    "post_diagnostics_us",
    "total_to_u_us",
    "total_function_us",
)


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--control-period-us", type=float, default=125.0)
    parser.add_argument(
        "--target-speed-factor",
        type=float,
        default=1.0,
        help="Target CPU single-thread speed divided by this PC's speed.",
    )
    parser.add_argument("--summary", type=Path, default=Path("u_timing_summary.csv"))
    parser.add_argument("--plot", type=Path, default=Path("u_timing.png"))
    args = parser.parse_args()
    if args.control_period_us <= 0.0 or args.target_speed_factor <= 0.0:
        parser.error("control period and target speed factor must be positive")

    with args.input_csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise RuntimeError("Timing CSV has no samples")
    missing = [name for name in PHASE_COLUMNS if name not in rows[0]]
    if missing:
        raise RuntimeError(f"Missing columns: {', '.join(missing)}")

    total_mean = statistics.fmean(float(row["total_function_us"]) for row in rows)
    summaries: list[dict[str, float | str]] = []
    for name in PHASE_COLUMNS:
        values = [float(row[name]) for row in rows]
        mean_us = statistics.fmean(values)
        summaries.append(
            {
                "phase": name,
                "mean_us": mean_us,
                "p50_us": percentile(values, 0.50),
                "p95_us": percentile(values, 0.95),
                "p99_us": percentile(values, 0.99),
                "max_us": max(values),
                "mean_share_of_function_pct": 100.0 * mean_us / total_mean,
                "target_cpu_p99_us_estimate": percentile(values, 0.99) / args.target_speed_factor,
            }
        )

    with args.summary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    total_to_u = [float(row["total_to_u_us"]) for row in rows]
    local_p99 = percentile(total_to_u, 0.99)
    target_p99 = local_p99 / args.target_speed_factor
    print(f"samples: {len(rows)}")
    print(f"local total-to-U p99: {local_p99:.3f} us")
    print(f"local deadline utilization: {100.0 * local_p99 / args.control_period_us:.2f}%")
    print(f"target CPU estimated p99: {target_p99:.3f} us")
    print(f"target deadline utilization: {100.0 * target_p99 / args.control_period_us:.2f}%")
    print(f"summary: {args.summary.resolve()}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is unavailable; skipped plot")
        return 0

    phase_rows = summaries[:8]
    labels = [str(row["phase"]).removesuffix("_us") for row in phase_rows]
    means = [float(row["mean_us"]) for row in phase_rows]
    p99s = [float(row["p99_us"]) for row in phase_rows]
    figure, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    axes[0].bar(labels, means, label="mean")
    axes[0].scatter(labels, p99s, color="black", s=20, label="p99", zorder=3)
    axes[0].set_ylabel("Time (us)")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].legend()
    axes[0].set_title("MPC U-computation stage latency")
    axes[1].plot(range(len(total_to_u)), total_to_u, linewidth=0.8)
    axes[1].axhline(args.control_period_us, color="red", linestyle="--", label="control period")
    axes[1].set_xlabel("Sample")
    axes[1].set_ylabel("Total to U (us)")
    axes[1].legend()
    figure.savefig(args.plot, dpi=180)
    print(f"plot: {args.plot.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
