#!/usr/bin/env bash
# 🧪📤 Build a CAM -usr_src directory: interface/cesm215 plus the column dump hooks.
# 🎛️🔕 interface/cesm215 itself stays untouched; the dump is only active with pumas_dump.nml in RUNDIR.
# Usage: make_src.sh <out_dir>
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${here}/../../.." && pwd)"
out="${1:?usage: make_src.sh <out_dir>}"
mkdir -p "${out}"
cp -L "${repo}"/interface/cesm215/*.F90 "${out}/"
cp "${repo}/tests/column/pumas_column_io.F90" "${out}/"
python3 "${here}/patch_micro_mg_cam.py" "${repo}/interface/cesm215/micro_mg_cam.F90" \
    "${out}/micro_mg_cam.F90" >/dev/null
echo "$(cd "${out}" && pwd)"
