"""Benchmark the instrumented DISCON DLL and report the MPC U-compute latency."""

from __future__ import annotations

import argparse
import csv
import ctypes
import math
import statistics
import time
from pathlib import Path


TIMING_NAMES = (
    "start_timestamp_us",
    "end_timestamp_us",
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
    parser.add_argument("--dll", type=Path, required=True)
    parser.add_argument("--in-file", type=Path, default=Path(__file__).resolve().parents[1] / "DISCON_MPC_CPP.IN")
    parser.add_argument("--warmup", type=int, default=20)
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--dt", type=float, default=0.001, help="DLL call interval in simulation seconds")
    parser.add_argument("--output", type=Path, default=Path("benchmark_u_timing_samples.csv"))
    args = parser.parse_args()
    if args.dt <= 0 or not math.isfinite(args.dt) or args.samples < 1 or args.warmup < 0:
        parser.error("dt and samples must be positive; warmup must be nonnegative")

    dll_path = args.dll.resolve()
    in_path = args.in_file.resolve()
    settings = {}
    for line in in_path.read_text(encoding="utf-8-sig").splitlines():
        line = line.split("#", 1)[0].split("!", 1)[0].strip()
        if line:
            key, value = line.split("=", 1)
            settings[key.strip()] = value.strip()
    fracture_time = float(settings["FRACTURE_TIME"])
    control_dt = float(settings["CTRL_DT"])
    max_calls = max(100, (args.warmup + args.samples) * math.ceil(control_dt / args.dt) * 10)
    controller = ctypes.CDLL(str(dll_path))
    discon = controller.DISCON
    discon.argtypes = [
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_char_p,
        ctypes.c_char_p,
        ctypes.c_char_p,
    ]
    get_timing = controller.MPC_GET_LAST_TIMING
    get_timing.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int]
    get_timing.restype = ctypes.c_int
    get_solve_info = controller.MPC_GET_LAST_SOLVE_INFO
    get_solve_info.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_int]
    get_solve_info.restype = ctypes.c_int

    swap = (ctypes.c_float * 3000)()
    fail = ctypes.c_int(0)
    message = ctypes.create_string_buffer(1024)
    in_bytes = str(in_path).encode("utf-8")
    out_bytes = str((Path.cwd() / "benchmark_u_timing.out").resolve()).encode("utf-8")
    swap[48] = len(message)
    swap[49] = len(in_bytes)
    swap[50] = len(out_bytes)
    swap[3] = 0.0
    swap[19] = 120.0
    swap[20] = 1.1
    swap[26] = 12.0
    swap[1018] = 0.0
    swap[1019] = 0.0
    swap[1020] = 0.0

    init_start = time.perf_counter_ns()
    swap[0] = 0.0
    swap[1] = 0.0
    discon(swap, ctypes.byref(fail), in_bytes, out_bytes, message)
    init_us = (time.perf_counter_ns() - init_start) / 1000.0
    if fail.value < 0:
        raise RuntimeError(message.value.decode("utf-8", errors="replace"))

    rows: list[list[float]] = []
    solve_info_rows: list[list[int]] = []
    call_index = 0
    solve_index = 0
    while solve_index < args.warmup + args.samples:
        call_index += 1
        if call_index > max_calls:
            detail = message.value.decode("utf-8", errors="replace")
            raise RuntimeError(f"No sufficient timing samples after {max_calls} calls. Enable MPC_ENABLE_TIMING and check solver failures. {detail}")
        swap[0] = 1.0
        swap[1] = fracture_time + call_index * args.dt
        discon(swap, ctypes.byref(fail), in_bytes, out_bytes, message)
        if fail.value < 0:
            raise RuntimeError(message.value.decode("utf-8", errors="replace"))
        timing = (ctypes.c_double * len(TIMING_NAMES))()
        count = get_timing(timing, len(TIMING_NAMES))
        if count == 0:
            continue
        if count != len(TIMING_NAMES):
            detail = message.value.decode("utf-8", errors="replace")
            raise RuntimeError(f"No successful timing sample at call {call_index}: {detail}")
        solve_index += 1
        if solve_index > args.warmup:
            rows.append(list(timing))
            solve_info = (ctypes.c_int * 2)()
            if get_solve_info(solve_info, 2) != 2:
                raise RuntimeError("Could not read QP solve mode")
            solve_info_rows.append(list(solve_info))

    # Exercise the normal Bladed/OpenFAST final-status path so the DLL flushes
    # its buffered per-call timing CSV as it would at the end of a simulation.
    swap[0] = -1.0
    swap[1] = fracture_time + (call_index + 1) * args.dt
    discon(swap, ctypes.byref(fail), in_bytes, out_bytes, message)

    print(f"DLL: {dll_path}")
    print(f"IN: {in_path}")
    print(f"one-time initialization: {init_us:.3f} us")
    print(f"samples: {len(rows)} after {args.warmup} warm-up solves")
    print(f"controller calls: {call_index} at dt={args.dt:.6f} s")
    output_path = args.output.resolve()
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("qp_solve_mode", "nwsr_used", *TIMING_NAMES))
        writer.writerows((*solve_info_rows[index], *row) for index, row in enumerate(rows))
    print(f"raw samples: {output_path}")
    mode_counts = {mode: sum(info[0] == mode for info in solve_info_rows) for mode in (0, 1, 2)}
    print(f"solve modes (cold/hot/hot-failed-cold): {mode_counts[0]}/{mode_counts[1]}/{mode_counts[2]}")
    print("phase,mean_us,p50_us,p95_us,p99_us,max_us")
    for column, name in enumerate(TIMING_NAMES[2:], start=2):
        values = [row[column] for row in rows]
        print(
            f"{name},{statistics.fmean(values):.3f},{percentile(values, 0.50):.3f},"
            f"{percentile(values, 0.95):.3f},{percentile(values, 0.99):.3f},{max(values):.3f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
