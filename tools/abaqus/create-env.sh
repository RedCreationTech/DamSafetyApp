#!/usr/bin/env bash

set -euo pipefail

readonly PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly ENV_PREFIX="${PROJECT_DIR}/.build/abaqus-converter-env"
readonly MAMBA_BIN="/home/kevin/miniforge3/bin/mamba"

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  "${MAMBA_BIN}" create -y -p "${ENV_PREFIX}" \
    'python=3.11' 'numpy=2.3' 'netcdf4=1.7'
fi

"${ENV_PREFIX}/bin/python" -c \
  'import netCDF4, numpy; print(f"netCDF4={netCDF4.__version__} numpy={numpy.__version__}")'
