#!/usr/bin/env bash

set -euo pipefail
source "$(dirname "$0")/common.sh"

mkdir -p "${PROJECT_DIR}/.upstream"

if [[ ! -d "${BLACKBEAR_DIR}/.git" ]]; then
  git clone --filter=blob:none --no-checkout \
    https://github.com/idaholab/blackbear.git "$BLACKBEAR_DIR"
fi

git -C "$BLACKBEAR_DIR" fetch --filter=blob:none origin "$BLACKBEAR_SHA"
git -C "$BLACKBEAR_DIR" checkout --detach "$BLACKBEAR_SHA"
git -C "$BLACKBEAR_DIR" submodule update --init contrib/neml moose

assert_checkout "$BLACKBEAR_DIR" "$BLACKBEAR_SHA"
assert_checkout "$MOOSE_DIR" "$MOOSE_SHA"
assert_checkout "$NEML_DIR" "$NEML_SHA"

printf 'BlackBear %s\nMOOSE %s\nNEML %s\n' \
  "$BLACKBEAR_SHA" "$MOOSE_SHA" "$NEML_SHA"
