"""Extract the six local 5 MW MPC aerodynamic Jacobian tables in SI units.

The controller model needs partial derivatives of aerodynamic torque and thrust
with respect to rotor speed, collective pitch, and wind speed while all other
OpenFAST states are held fixed. These are C/D Jacobian entries, not full-model
steady-state gains of the form -C A^-1 B + D.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_WINDS = tuple(range(2, 31))
DEFAULT_SPEEDS = tuple(0.5 * index for index in range(2, 25))
FILE_PATTERN = re.compile(r"(?P<wind>\d+(?:\.\d+)?)mps_(?P<rpm>\d+(?:\.\d+)?)rpm\.lin$")
TABLE_DEFINITIONS = {
    "GTomega.csv": ("ADRtAeroMxh_[Nm]", "C", "d_psi_rot_[rad/s]", "Nm/(rad/s)"),
    "GTbeta.csv": ("ADRtAeroMxh_[Nm]", "D", "PitchColl_[rad]", "Nm/rad"),
    "GTu.csv": ("ADRtAeroMxh_[Nm]", "D", "WS_[m/s]", "Nm/(m/s)"),
    "GFomega.csv": ("ADRtAeroFxh_[N]", "C", "d_psi_rot_[rad/s]", "N/(rad/s)"),
    "GFbeta.csv": ("ADRtAeroFxh_[N]", "D", "PitchColl_[rad]", "N/rad"),
    "GFu.csv": ("ADRtAeroFxh_[N]", "D", "WS_[m/s]", "N/(m/s)"),
}


def parse_values(raw: str) -> tuple[float, ...]:
    return tuple(float(item.strip()) for item in raw.split(",") if item.strip())


def format_point(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def find_unique_column(columns, preferred: str, contains: tuple[str, ...]) -> str:
    if preferred in columns:
        return preferred
    matches = [str(name) for name in columns if all(token.lower() in str(name).lower() for token in contains)]
    if len(matches) != 1:
        raise KeyError(f"Expected one column matching {contains}, found {matches}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lin-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, help="Manifest from generate_5mw_linearization_cases.py")
    parser.add_argument(
        "--override-lin-root",
        type=Path,
        action="append",
        default=[],
        help="Flat {wind}mps_{rpm}rpm.lin directory whose files replace matching base-grid points.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("generated_5mw_gain_tables"))
    parser.add_argument("--config-template", type=Path, help="Write a matching DISCON .IN beside the tables")
    parser.add_argument(
        "--toolbox-root",
        type=Path,
        default=Path(r"D:\zju_soft_study\cangnan_fast\openfast_toolbox"),
    )
    parser.add_argument("--wind-points")
    parser.add_argument("--speed-points")
    args = parser.parse_args()

    sys.path.insert(0, str(args.toolbox_root.resolve()))
    from openfast_toolbox.io.fast_linearization_file import FASTLinearizationFile

    manifest = pd.read_csv(args.manifest.resolve()) if args.manifest else None
    if manifest is not None:
        required = {"wind_speed_mps", "target_rpm", "expected_lin_file"}
        missing_manifest_columns = sorted(required - set(manifest.columns))
        if missing_manifest_columns:
            raise KeyError(f"Manifest missing columns: {missing_manifest_columns}")
    wind_points = parse_values(args.wind_points) if args.wind_points else (
        tuple(sorted(manifest["wind_speed_mps"].astype(float).unique())) if manifest is not None else DEFAULT_WINDS
    )
    speed_points = parse_values(args.speed_points) if args.speed_points else (
        tuple(sorted(manifest["target_rpm"].astype(float).unique())) if manifest is not None else DEFAULT_SPEEDS
    )
    wind_index = {value: index for index, value in enumerate(wind_points)}
    speed_index = {value: index for index, value in enumerate(speed_points)}
    tables = {
        name: np.full((len(wind_points), len(speed_points)), np.nan, dtype=float)
        for name in TABLE_DEFINITIONS
    }

    files: dict[tuple[float, float], Path] = {}
    if manifest is not None:
        for _, row in manifest.iterrows():
            key = (float(row["wind_speed_mps"]), float(row["target_rpm"]))
            path = Path(row["expected_lin_file"]).resolve()
            if key in files:
                raise RuntimeError(f"Duplicate manifest operating point {key}")
            files[key] = path
    else:
        for path in args.lin_root.resolve().glob("*.lin"):
            match = FILE_PATTERN.fullmatch(path.name)
            if not match:
                continue
            key = (float(match.group("wind")), float(match.group("rpm")))
            if key in files:
                raise RuntimeError(f"Duplicate operating point {key}: {files[key]} and {path}")
            files[key] = path
        for override_root in args.override_lin_root:
            for path in override_root.resolve().glob("*.lin"):
                match = FILE_PATTERN.fullmatch(path.name)
                if not match:
                    continue
                key = (float(match.group("wind")), float(match.group("rpm")))
                files[key] = path

    expected = {(wind, speed) for wind in wind_points for speed in speed_points}
    missing = sorted(expected - set(files))
    extras = sorted(set(files) - expected)
    if missing:
        preview = ", ".join(f"{wind:g}mps/{rpm:g}rpm" for wind, rpm in missing[:20])
        raise RuntimeError(f"Missing {len(missing)} required .lin points: {preview}")
    missing_files = [path for key, path in files.items() if key in expected and not path.is_file()]
    if missing_files:
        preview = ", ".join(str(path) for path in missing_files[:10])
        raise FileNotFoundError(f"Manifest references {len(missing_files)} missing .lin files: {preview}")
    if extras:
        print(f"Ignoring {len(extras)} .lin files outside the requested grid")

    audit_rows = []
    for count, (wind, rpm) in enumerate(sorted(expected), start=1):
        path = files[(wind, rpm)]
        linearization = FASTLinearizationFile(str(path))
        frames = linearization.toDataFrame()
        c_matrix = frames["C"]
        d_matrix = frames["D"]
        state_omega = find_unique_column(c_matrix.columns, "d_psi_rot_[rad/s]", ("psi_rot", "rad/s"))
        input_pitch = find_unique_column(d_matrix.columns, "PitchColl_[rad]", ("pitchcoll",))
        input_wind = find_unique_column(d_matrix.columns, "WS_[m/s]", ("ws_", "m/s"))

        resolved = {
            "GTomega.csv": ("ADRtAeroMxh_[Nm]", state_omega),
            "GTbeta.csv": ("ADRtAeroMxh_[Nm]", input_pitch),
            "GTu.csv": ("ADRtAeroMxh_[Nm]", input_wind),
            "GFomega.csv": ("ADRtAeroFxh_[N]", state_omega),
            "GFbeta.csv": ("ADRtAeroFxh_[N]", input_pitch),
            "GFu.csv": ("ADRtAeroFxh_[N]", input_wind),
        }
        values = {}
        for table_name, (output_name, coordinate_name) in resolved.items():
            frame = c_matrix if TABLE_DEFINITIONS[table_name][1] == "C" else d_matrix
            if output_name not in frame.index:
                raise KeyError(f"{output_name} missing from {path}")
            value = float(frame.loc[output_name, coordinate_name])
            if not np.isfinite(value):
                raise ValueError(f"Non-finite {table_name} value in {path}")
            tables[table_name][wind_index[wind], speed_index[rpm]] = value
            values[table_name.removesuffix(".csv")] = value

        y = frames.get("y")
        actual_wind = float("nan")
        actual_rpm = float("nan")
        if y is not None and not y.empty:
            if "Wind1VelX_[m/s]" in y.columns:
                actual_wind = float(y.iloc[0]["Wind1VelX_[m/s]"])
            if "RotSpeed_[rpm]" in y.columns:
                actual_rpm = float(y.iloc[0]["RotSpeed_[rpm]"])
        audit_rows.append(
            {
                "grid_wind_mps": wind,
                "grid_rotor_speed_rpm": rpm,
                "actual_wind_mps": actual_wind,
                "actual_rotor_speed_rpm": actual_rpm,
                "wind_error_mps": actual_wind - wind,
                "rotor_speed_error_rpm": actual_rpm - rpm,
                "lin_time_s": float(linearization["t"]),
                "lin_file": str(path),
                **values,
            }
        )
        if count % 25 == 0 or count == len(expected):
            print(f"Parsed {count}/{len(expected)}", flush=True)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        if not np.isfinite(table).all():
            raise RuntimeError(f"Incomplete table after extraction: {name}")
        pd.DataFrame(table).to_csv(output_dir / name, index=False, header=False, float_format="%.10e")

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(output_dir / "gain_point_audit.csv", index=False, encoding="utf-8-sig")
    summary_rows = []
    for name, (output_name, matrix_name, coordinate_name, unit) in TABLE_DEFINITIONS.items():
        summary_rows.append(
            {
                "table": name,
                "definition": f"d({output_name})/d({coordinate_name}) from {matrix_name}",
                "unit": unit,
                "minimum": float(np.min(tables[name])),
                "maximum": float(np.max(tables[name])),
            }
        )
    pd.DataFrame(summary_rows).to_csv(output_dir / "table_summary.csv", index=False)
    (output_dir / "DISCON_grid_snippet.txt").write_text(
        f"{len(speed_points)}\n{len(wind_points)}\n"
        + ",".join(format_point(v) for v in speed_points)
        + "\n"
        + ",".join(format_point(v) for v in wind_points)
        + "\n",
        encoding="ascii",
    )
    if args.config_template:
        config_lines = args.config_template.resolve().read_text(encoding="utf-8", errors="strict").splitlines(True)
        data_indices = [
            index
            for index, line in enumerate(config_lines)
            if line.strip() and not line.lstrip().startswith("!")
        ]
        if len(data_indices) < 10:
            raise RuntimeError("Controller config template has fewer than 10 data lines")
        replacements = (
            str(len(speed_points)),
            str(len(wind_points)),
            ",".join(format_point(v) for v in speed_points),
            ",".join(format_point(v) for v in wind_points),
            "GTomega.csv",
            "GTbeta.csv",
            "GTu.csv",
            "GFomega.csv",
            "GFbeta.csv",
            "GFu.csv",
        )
        for data_index, value in zip(data_indices[:10], replacements):
            newline = "\r\n" if config_lines[data_index].endswith("\r\n") else "\n"
            config_lines[data_index] = value + newline
        (output_dir / "DISCON_MPC_CPP.IN").write_text("".join(config_lines), encoding="utf-8")
    print(f"Saved six SI-unit tables and audit data to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
