#!/usr/bin/env bash
# ♻️🧪 Bit-for-bit test: short CESM 2.1.5 run (default 3 time steps) of a reference case, CAM restart compared with a stored baseline.
# 📘📂 Strategy and switch matrix: README.md and switches.csv in this directory.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
W="${W:-/cluster/work/climate/${USER}/work/CESM}"
CASE="${CASE:-cesm215_pumas.F2000climo.bfb_dev}"                 # reference case (CESM 3-like namelist, graupel on)
CASEDIR="${CASEDIR:-${W}/cesm215/cime/scripts/${CASE}}"
BUILD="${BUILD:-1}"                                              # BUILD=0 skips case.build (namelist-only change)
STOP_OPTION="${STOP_OPTION:-nsteps}"                             # ⏱️ run length; baselines are only valid for one length
STOP_N="${STOP_N:-3}"
BASELINES="${BASELINES:-/cluster/work/climate/${USER}/cesm_runs/opumas_bfb_baselines/${CASE}/${STOP_OPTION}${STOP_N}}"
WALLTIME="${WALLTIME:-00:30:00}"                                 # ⏳ short limit so the run jobs backfill quickly
CLONES="${CLONES:-/cluster/work/climate/${USER}/cesm_runs/opumas_bfb_clones}"   # 🐑 per-label --keepexe clones for par
PAR="${PAR:-15}"                                                 # 🚦 parallel run jobs (≤ 16 Slurm jobs including herdr-compute)
MARK='! >>> bfb.sh switch block (removed automatically) >>>'

usage() {
	cat >&2 <<USAGE
Usage: $0 <command> [label]      (label from tests/switches.csv, default: base)
  run    [label]   build, run STOP_N STOP_OPTION, wait, check       (one Slurm run job, a few min)
  par    <prefix>  run all labels starting with prefix ('all' = every row but base) in parallel, one --keepexe clone per label
                   (build the reference case first with 'run base'); logs in CLONES/<CASE>.<label>/bfb.log
  submit [label]   build and submit only; follow with 'check' once the job is gone
  check  [label]   compare the case's restart with the label's baseline (first time: store it)
  save   [label]   overwrite the label's baseline with the case's restart (after an intended change)
  list             labels, expectations and which baselines exist
Environment: CASE, CASEDIR, BASELINES, BUILD=0, STOP_OPTION=nsteps, STOP_N=3, WALLTIME, CLONES, PAR=15
USAGE
	exit 1
}

row() { awk -F, -v l="$1" 'NR>1 && $1==l' "${HERE}/switches.csv"; }
restart() { ls -t "$(cd "${CASEDIR}" && ./xmlquery --value RUNDIR)"/"${CASE}".cam.r.*.nc 2>/dev/null | head -1; }
cime_env() {
	module purge
	module load stack/.2024-04-silent
	module load python/3.9.18 perl/5.38.0 libiconv/1.17-zuuenhv libxml2/2.10.3-5uc2w2w subversion/1.14.2
}
# 🔎🎯 exact job-name match: clone names contain the reference case name
busy() { squeue -u "${USER}" -h -o '%j' | grep -qxF -e "run.${CASE}" -e "st_archive.${CASE}"; }
strip_switch() { sed -i "/^${MARK}\$/,\$d" "${CASEDIR}/user_nl_cam"; }
compare() { (source ~/conda && python3 "${HERE}/cmp_restart.py" "$1" "$2"); }

submit() {
	local label="$1" r nl jobs
	r="$(row "${label}")"
	[[ -n "${r}" ]] || { echo "ERROR: label '${label}' not in ${HERE}/switches.csv" >&2; exit 1; }
	if busy; then
		echo "ERROR: a job of ${CASE} is still queued or running" >&2; exit 1
	fi
	nl="$(cut -d, -f3- <<<"${r}")"
	cime_env
	cd "${CASEDIR}"
	strip_switch
	if [[ -n "${nl}" ]]; then
		{ echo "${MARK}"; tr ';' '\n' <<<"${nl//@REPO@/${REPO}}"; } >>user_nl_cam
	fi
	./xmlchange STOP_OPTION="${STOP_OPTION}",STOP_N="${STOP_N}",REST_OPTION="${STOP_OPTION}",REST_N="${STOP_N}"
	./xmlchange CONTINUE_RUN=FALSE,RESUBMIT=0,DOUT_S=FALSE
	./xmlchange JOB_WALLCLOCK_TIME="${WALLTIME}" --subgroup case.run
	if [[ "${BUILD}" == 1 ]]; then
		# 🧩🔁 CIME never regenerates these lists, so a new source file (micro_pumas_sppe.F90) would be missing
		rm -f "$(./xmlquery --value EXEROOT)"/atm/obj/{Srcfiles,Depends}
		sbatch --wait "${W}/build.proxy.sh"
		[[ "$(./xmlquery --value BUILD_COMPLETE)" == TRUE ]] || { echo "ERROR: build failed, see ${CASEDIR}/build.err" >&2; exit 1; }
	fi
	rm -f "$(./xmlquery --value RUNDIR)"/"${CASE}".cam.r.*.nc      # a stale restart must never pass as a new result
	echo "${label}" >"${CASEDIR}/.bfb_label"
	./case.submit
}

wait_job() {
	sleep 20
	while busy; do sleep 30; done
}

check() {
	local label="$1" new base expect ran
	ran="$(cat "${CASEDIR}/.bfb_label" 2>/dev/null || echo base)"
	[[ "${ran}" == "${label}" ]] || { echo "ERROR: the case last ran label '${ran}', not '${label}'" >&2; exit 1; }
	new="$(restart)"
	[[ -n "${new}" ]] || { echo "ERROR: no cam.r restart, the run failed: ${CASEDIR}/CaseStatus" >&2; exit 1; }
	mkdir -p "${BASELINES}"
	base="${BASELINES}/${label}.cam.r.nc"
	expect="$(row "${label}" | cut -d, -f2)"
	if [[ "${expect}" == differs && -f "${BASELINES}/base.cam.r.nc" ]]; then
		if compare "${new}" "${BASELINES}/base.cam.r.nc" >/dev/null; then
			echo "FAIL ${label}: identical to 'base', the switch has no effect"; exit 1
		fi
		echo "ok   ${label}: differs from 'base' (switch is wired)"
	elif [[ "${expect}" == same && "${label}" != base && -f "${BASELINES}/base.cam.r.nc" ]]; then
		# 🟰🔎 'same' rows (switch inactive in this configuration) must equal 'base'
		compare "${new}" "${BASELINES}/base.cam.r.nc" >/dev/null || { echo "FAIL ${label}: expected identical to 'base'"; exit 1; }
		echo "ok   ${label}: identical to 'base' (as expected)"
	fi
	if [[ -f "${base}" ]]; then
		if compare "${new}" "${base}"; then echo "PASS ${label}: bit-for-bit with ${base}"
		else echo "FAIL ${label}: not bit-for-bit with ${base}"; exit 1; fi
	else
		cp "${new}" "${base}"; echo "NEW  ${label}: baseline stored"; echo "${base}"
	fi
}

par() {
	# 🐑⚡ one --keepexe clone per label (shares the reference EXEROOT and BASELINES), PAR runs at a time
	local prefix="$1" l c srcroot
	[[ -f "${BASELINES}/base.cam.r.nc" ]] || { echo "ERROR: no base baseline, run '$0 run base' first" >&2; exit 1; }
	mkdir -p "${CLONES}"
	export STOP_OPTION STOP_N WALLTIME
	cime_env
	srcroot="$(cd "${CASEDIR}" && ./xmlquery --value SRCROOT)"
	for l in $(awk -F, -v p="${prefix}" 'NR>1 && (p=="all" ? $1!="base" : index($1,p)==1) {print $1}' "${HERE}/switches.csv"); do
		c="${CLONES}/${CASE}.${l}"
		[[ -d "${c}" ]] || "${CASEDIR}/../create_clone" --case "${c}" --clone "${CASEDIR}" --keepexe >/dev/null
		# 📂🔁 create_clone resets SRCROOT to CIMEROOT/.., losing the CAM worktree (its namelist_definition.xml)
		sed -i "s|id=\"SRCROOT\" value=\"[^\"]*\"|id=\"SRCROOT\" value=\"${srcroot}\"|" "${c}"/env_case.xml "${c}"/LockedFiles/env_case.xml
		echo "${l}"
	done | xargs -P "${PAR}" -I{} bash -c \
		'CASE="${0}.{}" CASEDIR="${1}/${0}.{}" BASELINES="${2}" BUILD=0 "${3}" run {} >"${1}/${0}.{}/bfb.log" 2>&1; grep -E "^(PASS|NEW|FAIL|ok|ERROR)" "${1}/${0}.{}/bfb.log" || echo "ERROR {}: see ${1}/${0}.{}/bfb.log"' \
		"${CASE}" "${CLONES}" "${BASELINES}" "${HERE}/bfb.sh"
}

cmd="${1:-}"; label="${2:-base}"
case "${cmd}" in
	run) trap strip_switch EXIT; submit "${label}"; wait_job; check "${label}" ;;
	submit) submit "${label}" ;;
	par) [[ -n "${2:-}" ]] || usage; par "$2" ;;
	check) check "${label}"; strip_switch ;;
	save) mkdir -p "${BASELINES}"; cp "$(restart)" "${BASELINES}/${label}.cam.r.nc"; echo "${BASELINES}/${label}.cam.r.nc" ;;
	list)
		while IFS=, read -r l e n; do
			[[ "${l}" == label ]] && continue
			printf '%-22s %-8s %-9s %s\n' "${l}" "${e}" "$([[ -f "${BASELINES}/${l}.cam.r.nc" ]] && echo baseline || echo '-')" "${n}"
		done <"${HERE}/switches.csv" ;;
	*) usage ;;
esac
