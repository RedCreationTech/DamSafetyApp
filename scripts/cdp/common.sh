#!/usr/bin/env bash

set -euo pipefail

readonly PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly MOOSE_DIR="${MOOSE_DIR:-/home/kevin/DamSafetyApp/.upstream/blackbear/moose}"
readonly BLACKBEAR_DIR="${BLACKBEAR_DIR:-/home/kevin/DamSafetyApp/.upstream/blackbear}"
readonly BUILD_ENV_PREFIX="${BUILD_ENV_PREFIX:-/home/kevin/DamSafetyApp/.build/env}"
readonly CONDA_BIN="${CONDA_BIN:-/home/kevin/miniforge3/bin/conda}"

readonly EXPECTED_MOOSE_SHA="4bce02d91b56c7ed845a5747df4d24f415592504"
readonly EXPECTED_BLACKBEAR_SHA="1c190fd3d2b5f06a3518923f550a0e0a90b015d4"

assert_checkout() {
  local directory="$1"
  local expected_sha="$2"
  local actual_sha

  actual_sha="$(git -C "$directory" rev-parse HEAD)"
  if [[ "$actual_sha" != "$expected_sha" ]]; then
    printf 'Unexpected commit in %s: expected %s, got %s\n' \
      "$directory" "$expected_sha" "$actual_sha" >&2
    exit 1
  fi
}

assert_build_environment() {
  [[ -x "$CONDA_BIN" ]] || {
    printf 'Missing conda executable: %s\n' "$CONDA_BIN" >&2
    exit 1
  }
  [[ -x "$BUILD_ENV_PREFIX/bin/python" ]] || {
    printf 'Missing build environment: %s\n' "$BUILD_ENV_PREFIX" >&2
    exit 1
  }
  assert_checkout "$MOOSE_DIR" "$EXPECTED_MOOSE_SHA"
  assert_checkout "$BLACKBEAR_DIR" "$EXPECTED_BLACKBEAR_SHA"
}

run_in_build_environment() {
  "$CONDA_BIN" run --no-capture-output -p "$BUILD_ENV_PREFIX" "$@"
}
