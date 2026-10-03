#!/usr/bin/env bash
# 🧮✅ Standalone compile check of the oPUMAS files and the CESM 2.1.5 overrides (about 1 minute).
# 🔁📦 Reuses the compile line and the .mod files of an already built case, so nothing is linked or run.
# 🖥️⛔ Run inside herdr-compute, never on a login node.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IFACE="${IFACE:-${HERE}/../interface/cesm215}"      # sources to check (set to a worktree's interface dir to check a branch)
CASE="${CASE:-cesm215_pumas.F2000climo.bfb_dev}"
EXEROOT="${EXEROOT:-/cluster/work/climate/${USER}/cesm_runs/${CASE}/bld}"
OBJ="${EXEROOT}/atm/obj"
# 🪜📐 Dependency order: PUMAS core first, then the CAM overrides.
FILES=(pumas_kinds pumas_gamma_function micro_pumas_sppe micro_pumas_utils micro_pumas_diags pumas_stochastic_collect_tau
	module_neural_net tau_neural_net_quantile ML_fixer_check micro_pumas_v1
	wv_sat_methods micro_dv_nml hetfrz_classnuc hetfrz_classnuc_cam ndrop microp_aero check_energy
	macrop_driver clubb_intr micro_mg_cam)

case "$(hostname)" in eu-login-*) echo "ERROR: login node, run inside herdr-compute" >&2; exit 1 ;; esac
IFACE="$(cd "${IFACE}" && pwd)"
# 🔎📜 Newest build log that compiled micro_mg_cam.F90 (an incremental build may have skipped it).
line=""
for log in $(ls -t "${EXEROOT}"/atm.bldlog.*.gz 2>/dev/null); do
	line="$(zcat "${log}" | grep mpif90 | grep -m1 '/micro_mg_cam.F90' || true)"
	[[ -n "${line}" ]] && break
done
[[ -n "${line}" ]] || { echo "ERROR: no micro_mg_cam.F90 compile line in ${EXEROOT}/atm.bldlog.*.gz (set CASE= or EXEROOT= to a built case)" >&2; exit 1; }
line="${line% *}"                                    # drop the source file

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 netcdf-c 2>/dev/null || true
B="$(mktemp -d "${TMPDIR:-/tmp}/opumas_compile.XXXXXX")"
trap 'rm -rf "${B}"' EXIT
# 🧭📂 Fresh .mod files go to B and win over the case's; relative -I paths are anchored at the case's obj dir.
line="$(sed -e "s# -I\.\./# -I${OBJ}/../#g" -e "s# -I\. # -I${OBJ} #g" \
	-e "s#^mpif90  *-c#mpif90 -c -I${B} -J${B} -fallow-argument-mismatch#" <<<"${line}")"

cd "${B}"
start=${SECONDS}
for f in "${FILES[@]}"; do
	if ${line} "${IFACE}/${f}.F90" -o "${B}/${f}.o" 2>"${B}/${f}.err"; then
		echo "  ok    ${f}"
	else
		echo "  FAIL  ${f}"; cat "${B}/${f}.err" >&2
		echo "COMPILE CHECK FAILED: ${IFACE}/${f}.F90" >&2; exit 1
	fi
done
echo "COMPILE CHECK OK: ${#FILES[@]} files from ${IFACE} in $((SECONDS - start)) s"
