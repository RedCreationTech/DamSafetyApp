#!/usr/bin/env bash

set -euo pipefail

readonly PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly ENV_PREFIX="${PROJECT_DIR}/.build/abaqus-converter-env"
readonly MAMBA_BIN="/home/kevin/miniforge3/bin/mamba"

if [[ ! -x "${ENV_PREFIX}/bin/python" ]]; then
  "${MAMBA_BIN}" create -y -p "${ENV_PREFIX}" \
    'python=3.11' 'numpy=2.3' 'netcdf4=1.7' 'scipy=1.17'
elif ! "${ENV_PREFIX}/bin/python" -c 'import scipy' 2>/dev/null; then
  "${MAMBA_BIN}" install -y -p "${ENV_PREFIX}" 'scipy=1.17'
fi

"${ENV_PREFIX}/bin/python" -c \
  'import netCDF4, numpy, scipy; print(f"netCDF4={netCDF4.__version__} numpy={numpy.__version__} scipy={scipy.__version__}")'
