#!/usr/bin/env bash
# 🧪📤 Submit the 1-day column-dump run; dumps go to <RUNDIR>/column_dump.
# 🕐🧊 Default: steps 24 and 47 (midnight and late evening UTC), every 50th chunk.
set -euo pipefail
W="${W:-/cluster/work/climate/${USER}/work/CESM}"
CASE="${CASE:-cesm215_pumas.F2000climo.column_dump}"
NSTEPS="${NSTEPS:-24,47}"
STRIDE="${STRIDE:-50}"
casedir="${W}/cesm215/cime/scripts/${CASE}"
source /etc/profile.d/lmod.sh 2>/dev/null || true
module purge
module load stack/.2024-04-silent
module load python/3.9.18 perl/5.38.0 libiconv/1.17-zuuenhv libxml2/2.10.3-5uc2w2w subversion/1.14.2
cd "${casedir}"
rundir="$(./xmlquery --value RUNDIR)"
mkdir -p "${rundir}/column_dump"
rm -f "${rundir}"/column_dump/*.nc
printf "&pumas_dump\n dump_dir = '%s'\n dump_nsteps = %s\n chunk_stride = %s\n/\n" \
	"${rundir}/column_dump" "${NSTEPS}" "${STRIDE}" >"${rundir}/pumas_dump.nml"
./case.submit
echo "${rundir}/column_dump"
