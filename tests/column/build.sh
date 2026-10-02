#!/usr/bin/env bash
# 🛠️📦 Build the standalone PUMAS column driver (about 30 s).
# 🧮⚖️ Default flags avoid FMA contraction and fast math so the Julia twin can match closely;
# OPT="-O2" reproduces the optimisation of the CESM 2.1.5 case build (x86-64 baseline, no FMA).
# 🖥️⛔ Run inside herdr-compute, never on a login node.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "${HERE}/../.." && pwd)"
IFACE="${SRC}/interface/cesm215"
BUILD="${BUILD:-${HERE}/build}"
OPT="${OPT:--O0}"
CAPTURE="${CAPTURE:-0}"   # 🔬 CAPTURE=1 builds pumas_column_capture, which records every process call
SPPE="${SPPE:-0}"         # 🧊 SPPE=1 builds pumas_column_sppe from the Fortran of branch ${SPPE_REF}
SPPE_REF="${SPPE_REF:-dev/Sppe}"
FILES=("${SRC}/pumas_kinds.F90" "${SRC}/pumas_gamma_function.F90" "${IFACE}/wv_sat_methods.F90"
	"${SRC}/micro_pumas_utils.F90" "${SRC}/micro_pumas_diags.F90"
	"${SRC}/pumas_stochastic_collect_tau.F90" "${SRC}/module_neural_net.F90"
	"${SRC}/tau_neural_net_quantile.F90" "${SRC}/ML_fixer_check.F90"
	"${SRC}/micro_pumas_v1.F90" "${HERE}/stubs.F90" "${HERE}/pumas_column_io.F90" "${HERE}/driver.F90")

case "$(hostname)" in eu-login-*) echo "ERROR: login node, run inside herdr-compute" >&2; exit 1 ;; esac
set +u
source /etc/profile.d/lmod.sh 2>/dev/null || true
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 netcdf-c netcdf-fortran
set -u

exe=pumas_column_driver
flags=(${OPT} -fno-fast-math -ffp-contract=off -ffree-line-length-none -g -fbacktrace $(nf-config --fflags))
if [[ "${CAPTURE}" == 1 ]]; then
	BUILD="${BUILD}/capture"
	exe=pumas_column_capture
	flags+=(-DCAPTURE)
	mkdir -p "${BUILD}"
	python3 "${HERE}/capture/patch_v1.py" "${SRC}" "${BUILD}" >/dev/null
	for i in "${!FILES[@]}"; do
		case "${FILES[$i]}" in
		"${SRC}/micro_pumas_utils.F90")
			FILES[$i]="${HERE}/capture/pumas_column_capture.F90 ${BUILD}/micro_pumas_utils.F90" ;;
		"${SRC}/micro_pumas_v1.F90") FILES[$i]="${BUILD}/micro_pumas_v1.F90" ;;
		esac
	done
	read -r -a FILES <<<"${FILES[*]}"
fi
if [[ "${SPPE}" == 1 ]]; then
	# 🧊🌿 dev/Sppe changes only these three files; take them read-only from the branch.
	BUILD="${BUILD}/sppe"
	exe=pumas_column_sppe
	flags+=(-DSPPE)
	mkdir -p "${BUILD}/src"
	git -C "${SRC}" archive "${SPPE_REF}" micro_pumas_sppe.F90 micro_pumas_utils.F90 \
		micro_pumas_v1.F90 | tar -x -C "${BUILD}/src"
	for i in "${!FILES[@]}"; do
		case "${FILES[$i]}" in
		"${SRC}/micro_pumas_utils.F90")
			FILES[$i]="${BUILD}/src/micro_pumas_sppe.F90 ${BUILD}/src/micro_pumas_utils.F90" ;;
		"${SRC}/micro_pumas_v1.F90") FILES[$i]="${BUILD}/src/micro_pumas_v1.F90" ;;
		esac
	done
	read -r -a FILES <<<"${FILES[*]}"
fi
mkdir -p "${BUILD}"
cd "${BUILD}"
for f in "${FILES[@]}"; do
	gfortran -c "${flags[@]}" "${f}" -o "$(basename "${f%.F90}").o"
done
gfortran -o "${exe}" *.o $(nf-config --flibs) -Wl,-rpath,"$(nf-config --prefix)/lib" \
	-Wl,-rpath,"$(nc-config --libdir)"
echo "${BUILD}/${exe}"
