#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly JOBS="${JOBS:-4}"
readonly METHOD="${METHOD:-opt}"

assert_build_environment
[[ -x "$PROJECT_DIR/DamSafetyApp-$METHOD" ]] || "$PROJECT_DIR/scripts/cdp/build-app.sh"
run_in_build_environment env MOOSE_DIR="$MOOSE_DIR" METHOD="$METHOD" \
  "$PROJECT_DIR/run_tests" -j "$JOBS" "$@"
