#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly TEST_JOBS="${TEST_JOBS:-2}"
readonly RUN_FULL="${RUN_FULL:-1}"
readonly TEST_LOG="${LOG_DIR}/TASK-MOOSE-004-official-gate.log"
readonly FULL_TEST_LOG="${LOG_DIR}/TASK-MOOSE-004-full-default.log"

assert_checkout "$BLACKBEAR_DIR" "$BLACKBEAR_SHA"
[[ -x "${BLACKBEAR_DIR}/blackbear-opt" ]] || {
  printf 'Missing blackbear-opt; run build-blackbear.sh first\n' >&2
  exit 1
}
mkdir -p "$LOG_DIR"

{
  printf 'BlackBear=%s\nMOOSE=%s\nNEML=%s\n' \
    "$BLACKBEAR_SHA" "$MOOSE_SHA" "$NEML_SHA"
  run_in_env "$BLACKBEAR_DIR/blackbear-opt" --version
  run_in_env "$BLACKBEAR_DIR/blackbear-opt" \
    --check-input \
    -i "$BLACKBEAR_DIR/test/tests/concrete_ASR_swelling/asr_confined.i"
  (
    cd "$BLACKBEAR_DIR"
    run_in_env env MOOSE_DIR="$MOOSE_DIR" ./run_tests \
      -j "$TEST_JOBS" \
      --re='concrete_ASR_(swelling|validation)'
  )
} 2>&1 | tee "$TEST_LOG"

if [[ "$RUN_FULL" == "1" ]]; then
  (
    cd "$BLACKBEAR_DIR"
    run_in_env env MOOSE_DIR="$MOOSE_DIR" ./run_tests -j "$TEST_JOBS"
  ) 2>&1 | tee "$FULL_TEST_LOG"
elif [[ "$RUN_FULL" != "0" ]]; then
  printf 'RUN_FULL must be 0 or 1\n' >&2
  exit 1
fi
