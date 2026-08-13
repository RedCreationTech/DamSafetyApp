#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly MPI_PROCS="${MPI_PROCS:-2}"
readonly AVAILABLE_CPUS="$(getconf _NPROCESSORS_ONLN)"
readonly INPUT_FILE="${BLACKBEAR_DIR}/test/tests/concrete_ASR_swelling/asr_confined.i"
readonly SMOKE_LOG="${LOG_DIR}/TASK-MOOSE-003-runtime-smoke.log"

[[ "$MPI_PROCS" =~ ^[1-9][0-9]*$ ]] || {
  printf 'MPI_PROCS must be a positive integer\n' >&2
  exit 1
}
(( MPI_PROCS <= AVAILABLE_CPUS )) || {
  printf 'Requested %s MPI ranks, but only %s logical CPUs are available\n' \
    "$MPI_PROCS" "$AVAILABLE_CPUS" >&2
  exit 1
}
[[ -x "${BLACKBEAR_DIR}/blackbear-opt" ]] || {
  printf 'Missing blackbear-opt; run build-blackbear.sh first\n' >&2
  exit 1
}
mkdir -p "$LOG_DIR"

{
  printf 'serial_check\n'
  run_in_env "$BLACKBEAR_DIR/blackbear-opt" --check-input -i "$INPUT_FILE"
  printf 'mpi_check ranks=%s\n' "$MPI_PROCS"
  run_in_env mpiexec -n "$MPI_PROCS" \
    "$BLACKBEAR_DIR/blackbear-opt" --check-input -i "$INPUT_FILE"
} 2>&1 | tee "$SMOKE_LOG"
