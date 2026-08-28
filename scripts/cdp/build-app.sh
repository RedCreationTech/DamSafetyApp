#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly JOBS="${JOBS:-4}"
readonly METHOD="${METHOD:-opt}"

assert_build_environment
run_in_build_environment env MOOSE_DIR="$MOOSE_DIR" METHOD="$METHOD" \
  make -C "$PROJECT_DIR" -j "$JOBS"

readonly EXECUTABLE="$PROJECT_DIR/DamSafetyApp-$METHOD"
[[ -x "$EXECUTABLE" ]] || {
  printf 'Build completed without executable: %s\n' "$EXECUTABLE" >&2
  exit 1
}

sha256sum "$EXECUTABLE"
