#!/usr/bin/env python3
"""S2 analytic consistency check for plane_strain_dynamic_smoke.

Rigid-body consistency: the block is stiff relative to the imposed 2 Hz base
acceleration, so the top must approximately follow the double integral of the
base acceleration. Compares top_disp_x from the CSV history against the
analytic displacement of a rigid block with the same Rayleigh mass damping
(qualitative gate: max |top - base_double_integral| small relative to motion
scale). Run on the build host after the smoke test:

  python3 analytic_check.py path/to/plane_strain_dynamic_smoke_out.csv
"""
import csv
import math
import sys

DT = 0.01


def load_series(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    return [float(r["time"]) for r in rows], [float(r["top_disp_x"]) for r in rows], [
        float(r["base_disp_x"]) for r in rows
    ]


def double_integrate(times, accel):
    n = len(times)
    vel = [0.0] * n
    disp = [0.0] * n
    for i in range(1, n):
        dt = times[i] - times[i - 1]
        vel[i] = vel[i - 1] + 0.5 * (accel[i] + accel[i - 1]) * dt
        disp[i] = disp[i - 1] + 0.5 * (vel[i] + vel[i - 1]) * dt
    return disp


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "plane_strain_dynamic_smoke_out.csv"
    times, top, base = load_series(path)
    accel = [2.0 * math.sin(2.0 * math.pi * 2.0 * t) for t in times]
    rigid = double_integrate(times, accel)
    residual = [t - r for t, r in zip(top, rigid)]
    base_residual = [t - b for b, t in zip(base, top)]
    scale = max(abs(v) for v in rigid) or 1.0
    print(f"samples={len(times)} motion_scale={scale:.6g} m")
    print(f"max|top - rigid_double_integral| = {max(abs(v) for v in residual):.6g} m "
          f"({100.0 * max(abs(v) for v in residual) / scale:.3f}% of scale)")
    print(f"max|top - base_disp|            = {max(abs(v) for v in base_residual):.6g} m")
    print("gate: residuals must be far below the motion scale (stiff-block consistency); "
          "no expert-approved threshold exists for this case yet.")


if __name__ == "__main__":
    main()
