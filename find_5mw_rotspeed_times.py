"""Find wind-specific post-fracture rotor-speed crossings for 5 MW linearization."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_TARGETS = tuple(0.5 * index for index in range(24, 1, -1))


def parse_values(raw: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in raw.split(",") if item.strip())


def find_column(frame: pd.DataFrame, base_name: str) -> str:
    matches = [name for name in frame.columns if str(name) == base_name or str(name).startswith(base_name + "_")]
    if len(matches) != 1:
        raise KeyError(f"Expected one {base_name} column, found {matches}")
    return str(matches[0])


def load_output(path: Path, toolbox_root: Path) -> pd.DataFrame:
    sys.path.insert(0, str(toolbox_root.resolve()))
    try:
        from openfast_toolbox.io.fast_output_file import FASTOutputFile
    except ImportError:
        from pyFAST.input_output import FASTOutputFile
    frame = FASTOutputFile(str(path)).toDataFrame()
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"Could not convert {path} to a dataframe")
    return frame


def smooth(values: np.ndarray, samples: int) -> np.ndarray:
    samples = max(1, samples | 1)
    if samples == 1:
        return values.copy()
    pad = samples // 2
    return np.convolve(np.pad(values, pad, mode="edge"), np.ones(samples) / samples, mode="valid")


def infer_wind(path: Path, frame: pd.DataFrame) -> float:
    match = re.search(r"(?P<wind>\d+(?:\.\d+)?)mps", str(path), re.IGNORECASE)
    if match:
        return float(match.group("wind"))
    wind_column = find_column(frame, "Wind1VelX")
    return float(np.nanmedian(frame[wind_column].to_numpy(dtype=float)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases-root", type=Path, required=True)
    parser.add_argument("--case-glob", default="steady_*mps")
    parser.add_argument("--out-glob", default="*.out")
    parser.add_argument("--fracture-time", type=float, default=30.0)
    parser.add_argument("--targets", default=",".join(f"{v:g}" for v in DEFAULT_TARGETS))
    parser.add_argument("--smooth-window-s", type=float, default=0.10)
    parser.add_argument("--local-window-s", type=float, default=0.30)
    parser.add_argument("--min-decel-rpmps", type=float, default=0.01)
    parser.add_argument("--points-output", type=Path, default=Path("5mw_linearization_candidates.csv"))
    parser.add_argument("--schedule-output", type=Path, default=Path("5mw_linearization_schedule.csv"))
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument(
        "--toolbox-root",
        type=Path,
        default=Path(r"D:\zju_soft_study\cangnan_fast\openfast_toolbox"),
    )
    args = parser.parse_args()

    targets = parse_values(args.targets)
    candidate_rows = []
    selected_rows = []
    missing_rows = []
    case_dirs = sorted(path for path in args.cases_root.resolve().glob(args.case_glob) if path.is_dir())
    if not case_dirs:
        raise FileNotFoundError("No time-domain case directories matched")

    for case_dir in case_dirs:
        outputs = sorted(case_dir.glob(args.out_glob))
        if len(outputs) != 1:
            raise RuntimeError(f"Expected one output in {case_dir}, found {len(outputs)}")
        out_file = outputs[0]
        frame = load_output(out_file, args.toolbox_root)
        time_col = find_column(frame, "Time")
        speed_col = find_column(frame, "RotSpeed")
        time_values = frame[time_col].to_numpy(dtype=float)
        speed_values = frame[speed_col].to_numpy(dtype=float)
        valid = np.isfinite(time_values) & np.isfinite(speed_values)
        time_values = time_values[valid]
        speed_values = speed_values[valid]
        order = np.argsort(time_values)
        time_values = time_values[order]
        speed_values = speed_values[order]
        dt = float(np.median(np.diff(time_values)))
        speed_smoothed = smooth(speed_values, int(round(args.smooth_window_s / dt)))
        speed_rate = np.gradient(speed_smoothed, time_values)
        speed_accel = np.gradient(speed_rate, time_values)
        wind = infer_wind(case_dir, frame)

        for target in targets:
            indices = np.where(
                (time_values[1:] > args.fracture_time)
                & (speed_smoothed[:-1] >= target)
                & (speed_smoothed[1:] < target)
                & (speed_rate[1:] <= -abs(args.min_decel_rpmps))
            )[0] + 1
            target_candidates = []
            for index in indices:
                y0, y1 = speed_smoothed[index - 1], speed_smoothed[index]
                alpha = 0.0 if abs(y1 - y0) < 1e-12 else (target - y0) / (y1 - y0)
                crossing_time = float(time_values[index - 1] + alpha * (time_values[index] - time_values[index - 1]))
                local = np.abs(time_values - crossing_time) <= 0.5 * args.local_window_s
                if np.count_nonzero(local) >= 3:
                    coeff = np.polyfit(time_values[local] - crossing_time, speed_smoothed[local], 1)
                    residual = speed_smoothed[local] - np.polyval(coeff, time_values[local] - crossing_time)
                    ripple = float(np.ptp(residual))
                else:
                    ripple = float("inf")
                row = {
                    "case_dir": str(case_dir),
                    "out_file": str(out_file),
                    "wind_speed_mps": wind,
                    "target_rpm": target,
                    "lin_time_s": crossing_time,
                    "actual_rpm": float(np.interp(crossing_time, time_values, speed_values)),
                    "smoothed_rpm": float(np.interp(crossing_time, time_values, speed_smoothed)),
                    "decel_rate_rpmps": float(speed_rate[index]),
                    "accel_rpmps2": float(speed_accel[index]),
                    "local_ripple_rpm": ripple,
                }
                row["score"] = ripple + 0.01 * abs(row["accel_rpmps2"])
                target_candidates.append(row)

            if not target_candidates:
                missing_rows.append((wind, target, case_dir))
                continue
            chosen = min(target_candidates, key=lambda row: (row["score"], row["lin_time_s"]))
            for row in target_candidates:
                row["selected"] = row is chosen
                candidate_rows.append(row)
            selected_rows.append(chosen.copy())

    if not selected_rows:
        for wind, target, case_dir in missing_rows:
            print(f"MISSING wind={wind:g} m/s target={target:g} rpm in {case_dir}")
        raise RuntimeError(
            "No post-fracture downward crossings were found. Check --fracture-time and the output trajectory."
        )

    schedule = pd.DataFrame(selected_rows).sort_values(["wind_speed_mps", "lin_time_s"])
    schedule["lin_order"] = schedule.groupby("wind_speed_mps").cumcount() + 1
    points = pd.DataFrame(candidate_rows).sort_values(["wind_speed_mps", "target_rpm", "lin_time_s"])
    points.to_csv(args.points_output.resolve(), index=False, encoding="utf-8-sig")
    schedule.to_csv(args.schedule_output.resolve(), index=False, encoding="utf-8-sig")

    for wind, target, case_dir in missing_rows:
        print(f"MISSING wind={wind:g} m/s target={target:g} rpm in {case_dir}")
    print(f"Selected {len(schedule)} points across {schedule.wind_speed_mps.nunique()} wind cases")
    if missing_rows and args.require_complete:
        raise RuntimeError(f"Missing {len(missing_rows)} wind/RPM crossings; no points were fabricated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
