#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly JOBS="${JOBS:-4}"
readonly BUILD_LOG="${LOG_DIR}/TASK-MOOSE-002-build-blackbear.log"

assert_checkout "$BLACKBEAR_DIR" "$BLACKBEAR_SHA"
assert_checkout "$MOOSE_DIR" "$MOOSE_SHA"
assert_checkout "$NEML_DIR" "$NEML_SHA"
mkdir -p "$LOG_DIR"

{
  printf '\n=== build start %s jobs=%s ===\n' "$(date --iso-8601=seconds)" "$JOBS"
  run_in_env env METHOD=opt make -C "$BLACKBEAR_DIR" -j "$JOBS"
  printf '=== build end %s ===\n' "$(date --iso-8601=seconds)"
} 2>&1 | tee -a "$BUILD_LOG"

[[ -x "${BLACKBEAR_DIR}/blackbear-opt" ]] || {
  printf 'Build completed without blackbear-opt\n' >&2
  exit 1
}

sha256sum "${BLACKBEAR_DIR}/blackbear-opt"
