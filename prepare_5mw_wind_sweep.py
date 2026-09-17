"""Create isolated 5 MW steady-wind fracture cases for the MPC DLL."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


WINDS = (12, 14, 16, 18, 20, 22, 24)
SKIP_SUFFIXES = {".out", ".outb", ".lin", ".png", ".log", ".yaml"}
SKIP_NAME_PARTS = (".sum", ".ech", ".chkp")


def ignore_generated(_directory: str, names: list[str]) -> set[str]:
    ignored = set()
    for name in names:
        path = Path(name)
        lower = name.lower()
        if path.suffix.lower() in SKIP_SUFFIXES or any(part in lower for part in SKIP_NAME_PARTS):
            ignored.add(name)
    return ignored


def replace_keyed_line(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^.*?\s+{re.escape(key)}\s+-.*$", re.MULTILINE)
    updated, count = pattern.subn(f"{value}   {key}       - configured by prepare_5mw_wind_sweep.py", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not uniquely replace {key}")
    return updated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-case", type=Path, required=True)
    parser.add_argument("--source-baseline", type=Path, required=True)
    parser.add_argument("--openfast", type=Path, required=True)
    parser.add_argument("--dll", type=Path, required=True)
    parser.add_argument("--in-file", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("openfast_5mw_wind_sweep"))
    args = parser.parse_args()

    output_root = args.output_root.resolve()
    baseline_out = output_root / "5MW_Baseline"
    output_root.mkdir(parents=True, exist_ok=True)
    shutil.copytree(args.source_baseline.resolve(), baseline_out, dirs_exist_ok=True, ignore=ignore_generated)
    shutil.copy2(args.openfast.resolve(), output_root / "OpenFAST5.0.0.exe")

    source_case = args.source_case.resolve()
    source_dir = source_case.parent
    source_servo = source_dir / "NRELOffshrBsline5MW_OC3Monopile_ServoDyn_fragment30_lin.dat"
    dll_path = args.dll.resolve().as_posix()
    in_path = args.in_file.resolve().as_posix()

    for wind in WINDS:
        case_dir = output_root / f"steady_{wind:02d}mps"
        case_dir.mkdir(parents=True, exist_ok=True)
        for source in source_dir.iterdir():
            if source.is_file() and source.suffix.lower() == ".dat":
                shutil.copy2(source, case_dir / source.name)

        fst_text = source_case.read_text(encoding="utf-8", errors="replace")
        fst_text = replace_keyed_line(
            fst_text,
            "InflowFile",
            f'"../5MW_Baseline/NRELOffshrBsline5MW_InflowWind_Steady{wind}mps.dat"',
        )
        fst_text = replace_keyed_line(fst_text, "Linearize", "False")
        fst_path = case_dir / f"mpc_fracture_{wind:02d}mps.fst"
        with fst_path.open("w", encoding="utf-8", newline="\r\n") as stream:
            stream.write(fst_text)

        servo_text = source_servo.read_text(encoding="utf-8", errors="replace")
        servo_text = replace_keyed_line(servo_text, "DLL_FileName", f'"{dll_path}"')
        servo_text = replace_keyed_line(servo_text, "DLL_InFile", f'"{in_path}"')
        with (case_dir / source_servo.name).open("w", encoding="utf-8", newline="\r\n") as stream:
            stream.write(servo_text)

    print(output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
