#!/usr/bin/env bash
# 🔤🧪 Build the column driver from the CF-CCPP translation of this working tree and check
#    that every switch row and the time loop are bit-for-bit with the PUMAS-name golden files.
# 📂 The translated copy lives in build/cf_ccpp (git-ignored); run run.sh first.
# 🖥️⛔ Run inside herdr-compute, never on a login node.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${here}/../.." && pwd)"
DATA="${DATA:-${here}/data}"
LANGUAGE="${LANGUAGE:-CF-CCPP}"
tmp="${here}/build/cf_ccpp"

rm -rf "${tmp}"
mkdir -p "${tmp}/tests" "${tmp}/interface"
cp "${REPO}"/*.F90 "${REPO}/KBARF_tau_kernel.dat" "${tmp}/"
cp -r "${REPO}/translation" "${tmp}/"
cp -r "${REPO}/interface/cesm215" "${tmp}/interface/"
mkdir -p "${tmp}/tests/column"
cp "${here}"/*.F90 "${here}"/*.sh "${tmp}/tests/column/"
python3 "${tmp}/translation/translate.py" --from PUMAS --to "${LANGUAGE}" --root "${tmp}" >/dev/null
BUILD="${tmp}/build" "${tmp}/tests/column/build.sh" >/dev/null
driver="${tmp}/build/pumas_column_driver"

source ~/conda >/dev/null 2>&1 || true
cd "${DATA}"
checked=0
while IFS=, read -r label _; do
	[[ "${label}" == label || ! -f "golden_${label}.nc" ]] && continue
	"${driver}" "settings_${label}.nml" golden_inputs.nc "${tmp}/${label}.nc" >/dev/null
	python3 "${here}/compare.py" --quiet "${tmp}/${label}.nc" "golden_${label}.nc" >/dev/null ||
		{ echo "❌ ${label} differs" >&2; exit 1; }
	checked=$((checked + 1))
done <"${here}/../switches.csv"
"${driver}" settings_base.nml golden_inputs.nc "${tmp}/loop.nc" 48 >/dev/null
python3 "${here}/compare.py" --quiet "${tmp}/loop.nc" golden_loop.nc >/dev/null ||
	{ echo "❌ loop differs" >&2; exit 1; }
echo "✅ ${LANGUAGE} driver bit-for-bit on ${checked} switch rows and the 48-step loop"
echo "${tmp}"
