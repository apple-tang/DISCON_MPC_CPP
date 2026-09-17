"""Generate isolated 5 MW Linearize=True cases from a selected crossing schedule."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd


GENERATED_SUFFIXES = {".out", ".outb", ".lin", ".log", ".png", ".yaml"}


def ignore_generated(_directory: str, names: list[str]) -> set[str]:
    return {
        name
        for name in names
        if Path(name).suffix.lower() in GENERATED_SUFFIXES or ".sum" in name.lower()
    }


def replace_keyed_line(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^.*?\s+{re.escape(key)}\s+-.*$", re.MULTILINE)
    updated, count = pattern.subn(f"{value}   {key}       - generated 5 MW linearization setting", text, count=1)
    if count != 1:
        raise RuntimeError(f"Could not uniquely replace {key}")
    return updated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--schedule", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=Path("openfast_5mw_linearization_cases"))
    parser.add_argument("--source-baseline", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    schedule = pd.read_csv(args.schedule.resolve())
    required = {"case_dir", "wind_speed_mps", "target_rpm", "lin_time_s", "lin_order"}
    missing = sorted(required - set(schedule.columns))
    if missing:
        raise KeyError(f"Schedule missing columns: {missing}")

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    if args.source_baseline:
        shutil.copytree(
            args.source_baseline.resolve(),
            output_root / "5MW_Baseline",
            dirs_exist_ok=True,
            ignore=ignore_generated,
        )

    manifest_rows = []
    for wind, group in schedule.groupby("wind_speed_mps", sort=True):
        source_dirs = {Path(value).resolve() for value in group["case_dir"]}
        if len(source_dirs) != 1:
            raise RuntimeError(f"Wind {wind:g} m/s maps to multiple source case directories")
        source_dir = source_dirs.pop()
        if "out_file" in group.columns:
            out_stems = {Path(value).stem for value in group["out_file"]}
            if len(out_stems) != 1:
                raise RuntimeError(f"Wind {wind:g} m/s maps to multiple output roots")
            fst_files = [source_dir / f"{out_stems.pop()}.fst"]
        else:
            fst_files = sorted(source_dir.glob("*.fst"))
        if len(fst_files) != 1 or not fst_files[0].is_file():
            raise RuntimeError(f"Could not identify one source .fst in {source_dir}: {fst_files}")

        target_dir = output_root / f"linearization_{wind:g}mps"
        if target_dir.exists() and not args.force:
            raise FileExistsError(f"{target_dir} exists; use --force only for generated cases")
        shutil.copytree(source_dir, target_dir, dirs_exist_ok=args.force, ignore=ignore_generated)

        ordered = group.sort_values("lin_time_s").reset_index(drop=True)
        times = ordered["lin_time_s"].astype(float).tolist()
        fst_text = fst_files[0].read_text(encoding="utf-8", errors="replace")
        fst_text = replace_keyed_line(fst_text, "Linearize", "True")
        fst_text = replace_keyed_line(fst_text, "NLinTimes", str(len(times)))
        fst_text = replace_keyed_line(fst_text, "LinTimes", ",".join(f"{value:.6f}" for value in times))
        fst_text = replace_keyed_line(fst_text, "TMax", f"{max(times) + 0.5:.6f}")
        target_fst = target_dir / f"linearization_{wind:g}mps.fst"
        with target_fst.open("w", encoding="utf-8", newline="\r\n") as stream:
            stream.write(fst_text)
        copied_original = target_dir / fst_files[0].name
        if copied_original.exists() and copied_original != target_fst:
            copied_original.unlink()

        for order, row in ordered.iterrows():
            manifest_rows.append(
                {
                    "case_dir": str(target_dir),
                    "fst_file": str(target_fst),
                    "lin_order": order + 1,
                    "lin_time_s": float(row["lin_time_s"]),
                    "wind_speed_mps": float(wind),
                    "target_rpm": float(row["target_rpm"]),
                    "expected_lin_file": str(target_dir / f"{target_fst.stem}.{order + 1}.lin"),
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = output_root / "linearization_manifest.csv"
    manifest.to_csv(manifest_path, index=False, encoding="utf-8-sig")
    print(f"Generated {manifest.case_dir.nunique()} cases and {len(manifest)} requested matrices")
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
