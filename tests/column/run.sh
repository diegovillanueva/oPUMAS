#!/usr/bin/env bash
# 🧪📦 Regenerate the column golden files from the selected CAM columns.
# 📘📂 See README.md in this directory.
#   run.sh select <dump_dir>   pick columns from a CAM dump into ${SOURCE} (once per dump)
#   run.sh                     build the driver, self-check against CAM, write golden_<label>.nc
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE="${SOURCE:-/cluster/work/climate/${USER}/opumas_column/source}"
DATA="${DATA:-${here}/data}"
SELF_RTOL="${SELF_RTOL:-1e-12}"
# 🎛️🚫 Rows that are not PUMAS switches at all.
SKIP="${SKIP:-wsub_default}"
REPO="$(cd "${here}/../.." && pwd)"

source ~/conda >/dev/null 2>&1 || true

if [[ "${1:-}" == select ]]; then
	dump="${2:?usage: run.sh select <dump_dir>}"
	python3 "${here}/select_columns.py" "${dump}" "${SOURCE}"
	python3 "${here}/synthetic.py" "${SOURCE}/inputs.nc" "${SOURCE}/synthetic_inputs.nc"
	# 🔤🧭 CAM renames rhmini_const to rhmini on use, but gfortran writes the original name.
	sed "s/^ *RHMINI_CONST=/ RHMINI=/" "${dump}/settings.nml" >"${SOURCE}/settings.nml"
	echo "${SOURCE}/settings.nml"
	exit 0
fi

"${here}/build.sh" >/dev/null
CAPTURE=1 "${here}/build.sh" >/dev/null
driver="${here}/build/pumas_column_driver"
capture="${here}/build/capture/pumas_column_capture"
mkdir -p "${DATA}"
cd "${DATA}"

# ♻️✅ Self-check: the driver reproduces CAM's own outputs for the dumped columns.
"${driver}" "${SOURCE}/settings.nml" "${SOURCE}/inputs.nc" "${DATA}/selfcheck.nc" >/dev/null
python3 "${here}/compare.py" --rtol "${SELF_RTOL}" "${DATA}/selfcheck.nc" "${SOURCE}/cam_outputs.nc" \
	| tee "${DATA}/selfcheck.log"

# 🧷📚 Golden inputs: CAM columns followed by synthetic ones.
python3 - "${SOURCE}/inputs.nc" "${SOURCE}/synthetic_inputs.nc" "${DATA}/golden_inputs.nc" <<'PY'
import sys
import xarray as xr
cam, syn = (xr.open_dataset(f) for f in sys.argv[1:3])
xr.concat([cam, syn[list(cam.data_vars)]], dim="col", data_vars="all").to_netcdf(sys.argv[3])
PY

# 🎛️📄 One golden file per in-scope switch row: settings.nml plus the row's overrides.
golden_rows() {
	local driver="$1" table="$2" skip_prefix="${3:-}" label expect nl
	while IFS=, read -r label expect nl; do
		[[ "${label}" == label ]] && continue
		[[ -n "${skip_prefix}" && "${label}" == "${skip_prefix}"* ]] && continue
		[[ " ${SKIP} " == *" ${label} "* ]] && { echo "skip ${label}"; continue; }
		sed '$d' "${SOURCE}/settings.nml" >"settings_${label}.nml"
		[[ -n "${nl}" ]] && tr ';' '\n' <<<"${nl//@REPO@/${REPO}}" >>"settings_${label}.nml"
		echo "/" >>"settings_${label}.nml"
		"${driver}" "settings_${label}.nml" golden_inputs.nc "golden_${label}.nc" >/dev/null
		echo "${DATA}/golden_${label}.nc"
	done <"${table}"
}
# 🧊🔀 The sppe_ rows of the CESM table also set CAM-side controls (hetfrz); the column has
#    its own Sppe table below.
golden_rows "${driver}" "${here}/../switches.csv" sppe_
# 🧊🎛️ Sppe rows: Sppe namelist values on top of settings.nml (column-only table).
golden_rows "${driver}" "${here}/switches_sppe.csv"
python3 "${here}/compare.py" --quiet golden_sppe_base.nc golden_base.nc
# ⏱️📈 Time loop: 432 steps of 600 s (3 days) of the base settings, forward Euler between
#    calls; the step count is stored as attribute nsteps for the Julia test.
nsteps="${LOOP_STEPS:-432}"
"${driver}" settings_base.nml golden_inputs.nc golden_loop.nc "${nsteps}" >/dev/null
python3 -c 'import sys, netCDF4; d = netCDF4.Dataset(sys.argv[1], "a"); d.nsteps = int(sys.argv[2]); d.close()' \
	golden_loop.nc "${nsteps}"
echo "${DATA}/golden_loop.nc"
# 🔬🧷 Process-level data: arguments of every process call in the base step.
PUMAS_CAPTURE_FILE="${DATA}/golden_processes.nc" "${capture}" "${SOURCE}/settings.nml" \
	golden_inputs.nc capture_base.nc >/dev/null
python3 "${here}/compare.py" --quiet capture_base.nc golden_base.nc   # capture must not change results
# 🌧️🧷 The Seifert-Beheng warm rain routines only run with micro_pumas_warm_rain='sb2001'.
PUMAS_CAPTURE_FILE="${DATA}/golden_processes_sb2001.nc" "${capture}" settings_sb2001.nml \
	golden_inputs.nc capture_sb2001.nc >/dev/null
python3 "${here}/compare.py" --quiet capture_sb2001.nc golden_sb2001.nc
echo "${DATA}/golden_processes.nc"
echo "${DATA}/golden_processes_sb2001.nc"
echo "${DATA}/golden_inputs.nc"
echo "${DATA}/selfcheck.log"
