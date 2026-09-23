#!/usr/bin/env python3
"""Verify that a MOOSE checkpoint restart preserves Abaqus-CDP constitutive state.

Runs one continuous loading history and the same history split into
"phase 1 + checkpoint -> phase 2 restart", then compares every postprocessor.
The split run must reproduce the continuous run step by step, which is only true if the
restored solution *and* the stateful material history (DamageT/C, kappa, plastic strain,
total strain) come back consistently.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TEST_DIR = REPO / "test/tests/abaqus_cdp_stress_update"
MATERIAL_DIR = REPO / "test/tests/cdp_material_table/data"
INPUTS = ("single_hex8_common.i", "restart_continuous.i", "restart_phase1.i", "restart_phase2.i")
MARKER = "__RESTART_BASE__"
IGNORED_COLUMNS = {"time"}
RELATIVE_TOLERANCE = 1e-8
ABSOLUTE_TOLERANCE = 1e-12


def run(binary: Path, working: Path, case: str, *extra: str) -> None:
    result = subprocess.run([str(binary), "-i", case, *extra], cwd=working,
                            capture_output=True, text=True)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-25:])
        raise RuntimeError(f"{case} exited {result.returncode}:\n{tail}")


def read_csv(working: Path, base: str) -> dict[float, dict[str, float]]:
    path = working / f"{base}.csv"
    if not path.is_file():
        raise RuntimeError(f"missing output {path.name}")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = {float(row["time"]): {key: float(value) for key, value in row.items()
                                     if key not in IGNORED_COLUMNS}
                for row in csv.DictReader(stream)}
    return {time: values for time, values in sorted(rows.items())}


def compare(reference: dict[float, dict[str, float]], candidate: dict[float, dict[str, float]],
            label: str) -> dict[str, object]:
    shared = sorted(set(reference) & set(candidate))
    if not shared:
        raise RuntimeError(f"{label}: no common time steps")
    columns = sorted(set(reference[shared[0]]) - IGNORED_COLUMNS)
    worst: dict[str, float] = {}
    failures: list[str] = []
    for column in columns:
        maximum = max(abs(candidate[time][column] - reference[time][column]) for time in shared)
        scale = max(abs(reference[time][column]) for time in shared)
        relative = maximum / scale if scale > 0 else (maximum if maximum > 0 else 0.0)
        worst[column] = relative
        if maximum > ABSOLUTE_TOLERANCE + RELATIVE_TOLERANCE * scale:
            failures.append(f"{column}: max abs difference {maximum:.3e} (scale {scale:.3e})")
    return {"label": label, "shared_steps": len(shared), "candidate_steps": len(candidate),
            "max_relative_difference": worst, "failures": failures,
            "passed": not failures and len(candidate) >= len(reference)}


def find_restart_base(working: Path, file_base: str) -> Path:
    candidates = sorted(working.glob(f"{file_base}*cp*/**/*-restart-*.rd")) + \
        sorted(working.glob(f"{file_base}*cp*/*-restart-*.rd")) + \
        sorted(working.glob(f"{file_base}*cp*/*.rd"))
    if not candidates:
        raise RuntimeError(f"no checkpoint data under {working} for file_base={file_base}; "
                           f"found directories: {[p.name for p in working.iterdir() if p.is_dir()]}")
    match = re.match(r"^(.+?)-restart", candidates[0].name) or re.match(r"^(\d+)", candidates[0].name)
    if not match:
        raise RuntimeError(f"cannot derive restart base from {candidates[0].name}")
    step = match.group(1)
    return (candidates[0].parent / step).resolve()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path, help="DamSafetyApp-opt (or -devel)")
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--keep", action="store_true", help="keep the scratch layout for inspection")
    args = parser.parse_args()

    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"solver binary not found: {binary}")
    working = args.work_dir.resolve()
    if working.exists():
        shutil.rmtree(working)
    (working / "abaqus_cdp_stress_update").mkdir(parents=True)
    shutil.copytree(MATERIAL_DIR, working / "cdp_material_table" / "data")
    for name in INPUTS:
        shutil.copy2(TEST_DIR / name, working / "abaqus_cdp_stress_update" / name)
    scratch = working / "abaqus_cdp_stress_update"

    run(binary, scratch, "restart_continuous.i")
    run(binary, scratch, "restart_phase1.i", "--Executioner/end_time", "0.5")
    base = find_restart_base(scratch, "restart_phase1")
    template = (scratch / "restart_phase2.i").read_text()
    if MARKER not in template:
        raise SystemExit(f"{ 'restart_phase2.i' } lost its {MARKER} placeholder")
    (scratch / "restart_phase2_run.i").write_text(template.replace(MARKER, str(base)))
    run(binary, scratch, "restart_phase2_run.i")

    continuous = read_csv(scratch, "restart_continuous_out")
    phase1 = read_csv(scratch, "restart_phase1_out")
    phase2 = read_csv(scratch, "restart_phase2_out")

    checks = [compare(continuous, {t: v for t, v in phase1.items() if t <= 0.5 + 1e-9},
                      "phase1 versus continuous (pre-restart)"),
              compare(continuous, {t: v for t, v in phase2.items() if t >= 0.5 - 1e-9},
                      "phase2 versus continuous (post-restart)")]
    final = continuous[max(continuous)]
    stateful = {key: final[key] for key in final if re.search(r"damage|kappa|stiffness", key)}
    exercised = any(value > 0 for key, value in stateful.items() if "damage" in key or "kappa" in key)
    summary = {
        "restart_base": str(base),
        "continuous_steps": len(continuous), "phase2_steps": len(phase2),
        "checks": checks,
        "final_state_of_reference_run": {key: value for key, value in sorted(stateful.items())},
        "constitutive_history_actually_exercised": exercised,
        "passed": all(check["passed"] for check in checks) and exercised and math.isfinite(
            max(value for values in continuous.values() for value in values.values())),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not args.keep:
        shutil.rmtree(working, ignore_errors=True)
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
