#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

readonly MAMBA_BIN="/home/kevin/miniforge3/bin/mamba"
assert_command "$MAMBA_BIN"

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  mkdir -p "$(dirname "$ENV_PREFIX")"
  "$MAMBA_BIN" create -y -p "$ENV_PREFIX" 'moose-dev=2026.07.30=mpich'
fi

run_in_env python --version
run_in_env mpiexec --version
run_in_env conda list --explicit
