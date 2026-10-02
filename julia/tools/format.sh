#!/usr/bin/env bash
# 🎨📏 Format src/ and test/ with JuliaFormatter (.JuliaFormatter.toml), or check them.
# 🧪 `tools/format.sh --check` changes nothing and fails if a file needs formatting or a
#    line is longer than 92 columns.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${here}"
"${here}/tools/julia.sh" --project="${here}/tools/format" -e 'using Pkg; Pkg.instantiate()' \
    >/dev/null
overwrite=true
[[ "${1:-}" == "--check" ]] && overwrite=false
"${here}/tools/julia.sh" --project="${here}/tools/format" -e "
    using JuliaFormatter
    clean = format([\"src\", \"test\"]; overwrite = ${overwrite})
    clean || println(\"🎨 files needed formatting\")
    exit(clean || ${overwrite} ? 0 : 1)"
long="$(awk 'length > 92 { print FILENAME ":" FNR }' src/*.jl src/*/*.jl test/*.jl)"
if [[ -n "${long}" ]]; then
    echo "📏 lines longer than 92 columns:"
    echo "${long}"
    exit 1
fi
echo "🎨✅ formatted, no line over 92 columns"
