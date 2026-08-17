#!/usr/bin/env bash

# TASK-MOOSE-007-A: 复现官方 ASR 受约束圆柱回归 test_strip
# 基线: BlackBear 锁定提交 concrete_ASR_swelling/asr_confined.i（Exodiff + asr_confined.cmp）

set -euo pipefail
source "$(dirname "$0")/../p0/common.sh"

readonly TEST_JOBS="${TEST_JOBS:-2}"
readonly TEST_LOG="${LOG_DIR}/TASK-MOOSE-007-A-asr-confined-strip.log"

assert_checkout "$BLACKBEAR_DIR" "$BLACKBEAR_SHA"
assert_checkout "$MOOSE_DIR" "$MOOSE_SHA"
assert_checkout "$NEML_DIR" "$NEML_SHA"
[[ -x "${BLACKBEAR_DIR}/blackbear-opt" ]] || {
  printf 'Missing blackbear-opt; run build-blackbear.sh first\n' >&2
  exit 1
}
mkdir -p "$LOG_DIR"

{
  printf 'DamSafetyApp=%s\nBlackBear=%s\nMOOSE=%s\nNEML=%s\n' \
    "$(git -C "$PROJECT_DIR" rev-parse HEAD)" \
    "$BLACKBEAR_SHA" "$MOOSE_SHA" "$NEML_SHA"
  run_in_env "$BLACKBEAR_DIR/blackbear-opt" --version
  (
    cd "$BLACKBEAR_DIR"
    run_in_env env MOOSE_DIR="$MOOSE_DIR" ./run_tests \
      -j "$TEST_JOBS" \
      --re="concrete_ASR_swelling\.ASR_swelling/test_strip\$"
  )
} 2>&1 | tee "$TEST_LOG"
