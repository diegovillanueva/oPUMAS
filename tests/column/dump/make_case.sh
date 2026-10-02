#!/usr/bin/env bash
# 🧪📤 Create (once) and build the CAM column-dump case: a clone of bfb_dev with the dump hooks.
# 📘📂 Then submit with run_case.sh; outputs land in <RUNDIR>/column_dump.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${W:-/cluster/work/climate/${USER}/work/CESM}"
SRC_CASE="${SRC_CASE:-cesm215_pumas.F2000climo.bfb_dev}"
CASE="${CASE:-cesm215_pumas.F2000climo.column_dump}"
scripts="${W}/cesm215/cime/scripts"
casedir="${scripts}/${CASE}"
srcroot="${W}/cesm215_pumas_worktrees/cesm215_pumas"

source /etc/profile.d/lmod.sh 2>/dev/null || true
module purge
module load stack/.2024-04-silent
module load python/3.9.18 perl/5.38.0 libiconv/1.17-zuuenhv libxml2/2.10.3-5uc2w2w subversion/1.14.2

if [[ ! -d "${casedir}" ]]; then
	"${W}/cesm215/cime/scripts/create_clone" --case "${casedir}" --clone "${scripts}/${SRC_CASE}"
	cd "${casedir}"
	# 🧭📂 create_clone resets SRCROOT; restore it in both env_case copies.
	sed -i "s#id=\"SRCROOT\" value=\"[^\"]*\"#id=\"SRCROOT\" value=\"${srcroot}\"#" \
		env_case.xml LockedFiles/env_case.xml
fi
cd "${casedir}"
"${here}/make_src.sh" "${casedir}/column_dump_src"
opts="$(./xmlquery --value CAM_CONFIG_OPTS)"
opts="$(sed -E "s#-usr_src [^ ]+#-usr_src ${casedir}/column_dump_src#" <<<"${opts}")"
./xmlchange CAM_CONFIG_OPTS="${opts}"
./xmlchange EXEROOT="${W%/work/CESM}/cesm_runs/${CASE}/bld",RUNDIR="${W%/work/CESM}/cesm_runs/${CASE}/run"
./xmlchange STOP_OPTION=ndays,STOP_N=1,CONTINUE_RUN=FALSE,RESUBMIT=0,DOUT_S=FALSE
./case.setup --reset >/dev/null
sbatch --wait "${W}/build.proxy.sh"
[[ "$(./xmlquery --value BUILD_COMPLETE)" == TRUE ]] || { echo "ERROR: build failed: ${casedir}/build.err" >&2; exit 1; }
echo "${casedir}"
