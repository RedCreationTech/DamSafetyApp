#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly JOBS="${JOBS:-4}"
readonly METHOD="${METHOD:-opt}"

assert_build_environment
run_in_build_environment env MOOSE_DIR="$MOOSE_DIR" METHOD="$METHOD" \
  make -C "$PROJECT_DIR/unit" -j "$JOBS"
run_in_build_environment env MOOSE_DIR="$MOOSE_DIR" METHOD="$METHOD" \
  "$PROJECT_DIR/unit/run_tests" "$@"
