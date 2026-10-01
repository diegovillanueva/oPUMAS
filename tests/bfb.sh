#!/usr/bin/env bash
# ♻️🧪 Bit-for-bit test: 1-day CESM 2.1.5 run of a reference case, CAM restart compared with a stored baseline.
# 📘📂 Strategy and switch matrix: README.md and switches.csv in this directory.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
W="${W:-/cluster/work/climate/${USER}/work/CESM}"
CASE="${CASE:-cesm215_pumas.F2000climo.bfb_dev}"                 # reference case (CESM 3-like namelist, graupel on)
CASEDIR="${CASEDIR:-${W}/cesm215/cime/scripts/${CASE}}"
BASELINES="${BASELINES:-/cluster/work/climate/${USER}/cesm_runs/opumas_bfb_baselines/${CASE}}"
BUILD="${BUILD:-1}"                                              # BUILD=0 skips case.build (namelist-only change)
MARK='! >>> bfb.sh switch block (removed automatically) >>>'

usage() {
	cat >&2 <<USAGE
Usage: $0 <command> [label]      (label from tests/switches.csv, default: base)
  run    [label]   build, run 1 day, wait, check                    (one Slurm run job, about 10 min)
  submit [label]   build and submit only; follow with 'check' once the job is gone
  check  [label]   compare the case's restart with the label's baseline (first time: store it)
  save   [label]   overwrite the label's baseline with the case's restart (after an intended change)
  list             labels, expectations and which baselines exist
Environment: CASE, CASEDIR, BASELINES, BUILD=0
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
strip_switch() { sed -i "/^${MARK}\$/,\$d" "${CASEDIR}/user_nl_cam"; }
compare() { (source ~/conda && python3 "${HERE}/cmp_restart.py" "$1" "$2"); }

submit() {
	local label="$1" r nl jobs
	r="$(row "${label}")"
	[[ -n "${r}" ]] || { echo "ERROR: label '${label}' not in ${HERE}/switches.csv" >&2; exit 1; }
	if squeue -u "${USER}" -h -o '%j' | grep -q "${CASE}"; then
		echo "ERROR: a job of ${CASE} is still queued or running" >&2; exit 1
	fi
	nl="$(cut -d, -f3- <<<"${r}")"
	cime_env
	cd "${CASEDIR}"
	strip_switch
	if [[ -n "${nl}" ]]; then
		{ echo "${MARK}"; tr ';' '\n' <<<"${nl//@REPO@/${REPO}}"; } >>user_nl_cam
	fi
	./xmlchange STOP_OPTION=ndays,STOP_N=1,CONTINUE_RUN=FALSE,RESUBMIT=0,DOUT_S=FALSE
	if [[ "${BUILD}" == 1 ]]; then
		sbatch --wait "${W}/build.proxy.sh"
		[[ "$(./xmlquery --value BUILD_COMPLETE)" == TRUE ]] || { echo "ERROR: build failed, see ${CASEDIR}/build.err" >&2; exit 1; }
	fi
	rm -f "$(./xmlquery --value RUNDIR)"/"${CASE}".cam.r.*.nc      # a stale restart must never pass as a new result
	echo "${label}" >"${CASEDIR}/.bfb_label"
	./case.submit
}

wait_job() {
	sleep 20
	while squeue -u "${USER}" -h -o '%j' | grep -q "${CASE}"; do sleep 30; done
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
	fi
	if [[ -f "${base}" ]]; then
		if compare "${new}" "${base}"; then echo "PASS ${label}: bit-for-bit with ${base}"
		else echo "FAIL ${label}: not bit-for-bit with ${base}"; exit 1; fi
	else
		cp "${new}" "${base}"; echo "NEW  ${label}: baseline stored"; echo "${base}"
	fi
}

cmd="${1:-}"; label="${2:-base}"
case "${cmd}" in
	run) trap strip_switch EXIT; submit "${label}"; wait_job; check "${label}" ;;
	submit) submit "${label}" ;;
	check) check "${label}"; strip_switch ;;
	save) mkdir -p "${BASELINES}"; cp "$(restart)" "${BASELINES}/${label}.cam.r.nc"; echo "${BASELINES}/${label}.cam.r.nc" ;;
	list)
		while IFS=, read -r l e n; do
			[[ "${l}" == label ]] && continue
			printf '%-22s %-8s %-9s %s\n' "${l}" "${e}" "$([[ -f "${BASELINES}/${l}.cam.r.nc" ]] && echo baseline || echo '-')" "${n}"
		done <"${HERE}/switches.csv" ;;
	*) usage ;;
esac
