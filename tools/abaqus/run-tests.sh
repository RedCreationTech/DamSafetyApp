#!/usr/bin/env bash

set -euo pipefail

readonly PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
readonly PYTHON="${PROJECT_DIR}/.build/abaqus-converter-env/bin/python"

[[ -x "$PYTHON" ]] || {
  printf 'Missing converter environment; run tools/abaqus/create-env.sh\n' >&2
  exit 1
}

cd "$PROJECT_DIR"
"$PYTHON" -m unittest discover -s test/tools -p 'test_*.py' -v
