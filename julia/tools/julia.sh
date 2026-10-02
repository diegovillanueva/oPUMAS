#!/usr/bin/env bash
# 🟣🚀 Run Julia 1.10 on the oPUMAS project with a clean library path.
# 🧩🔧 Euler's default LD_LIBRARY_PATH holds a gcc-8 libgfortran that breaks NetCDF_jll.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec env -u LD_LIBRARY_PATH /cluster/software/manual/julia/1.10.3/x86_64/bin/julia \
    --project="${here}" "$@"
