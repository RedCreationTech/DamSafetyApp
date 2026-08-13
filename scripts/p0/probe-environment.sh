#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly PROBE_LOG="${LOG_DIR}/TASK-MOOSE-003-environment.log"
mkdir -p "$LOG_DIR"

{
  printf 'timestamp=%s\n' "$(date --iso-8601=seconds)"
  printf 'host=%s\n' "$(hostname)"
  printf 'kernel=%s\n' "$(uname -srmo)"
  printf 'logical_cpus=%s\n' "$(getconf _NPROCESSORS_ONLN)"
  awk '/MemTotal|MemAvailable/ {print}' /proc/meminfo
  df -h "$PROJECT_DIR"
  run_in_env python --version
  run_in_env mpicc --version | head -1
  run_in_env mpicxx --version | head -1
  run_in_env mpiexec --version | head -4
  run_in_env conda list | grep -E '^moose-(build|dev|libmesh|petsc|tools|wasp)[[:space:]]'
} | tee "$PROBE_LOG"
