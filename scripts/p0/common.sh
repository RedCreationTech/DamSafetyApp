#!/usr/bin/env bash

set -euo pipefail

readonly PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly BLACKBEAR_DIR="${PROJECT_DIR}/.upstream/blackbear"
readonly MOOSE_DIR="${BLACKBEAR_DIR}/moose"
readonly NEML_DIR="${BLACKBEAR_DIR}/contrib/neml"
readonly ENV_PREFIX="${PROJECT_DIR}/.build/env"
readonly LOG_DIR="${PROJECT_DIR}/.build/logs"
readonly CONDA_BIN="/home/kevin/miniforge3/bin/conda"

readonly BLACKBEAR_SHA="1c190fd3d2b5f06a3518923f550a0e0a90b015d4"
readonly MOOSE_SHA="4bce02d91b56c7ed845a5747df4d24f415592504"
readonly NEML_SHA="a01a27b524a737b6746a840150f5acc2bace778e"

assert_command() {
  command -v "$1" >/dev/null 2>&1 || {
    printf 'Required command is missing: %s\n' "$1" >&2
    exit 1
  }
}

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

run_in_env() {
  [[ -x "${ENV_PREFIX}/bin/python" ]] || {
    printf 'Missing isolated environment: %s\n' "$ENV_PREFIX" >&2
    exit 1
  }
  "$CONDA_BIN" run --no-capture-output -p "$ENV_PREFIX" "$@"
}
