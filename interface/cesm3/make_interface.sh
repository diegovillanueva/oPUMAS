#!/usr/bin/env bash
# 🌍🔌 Sppe CAM-side patches for CESM 3 / CAM7 (written against CAM tag cam6_4_204).
# 📘📂 See README.md in this directory; these patches are checked, not compiled or run.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAG="${TAG:-cam6_4_204}"                                       # CAM tag the patches were written against
RAW="${RAW:-https://raw.githubusercontent.com/ESCOMP/CAM/${TAG}}"
FILES=(src/physics/cam/microp_aero.F90 src/physics/cam/hetfrz_classnuc.F90
	src/physics/cam/hetfrz_classnuc_cam.F90 src/physics/cam/clubb_intr.F90
	src/physics/cam/macrop_driver.F90 src/physics/cam/micro_pumas_cam.F90
	bld/namelist_files/namelist_definition.xml)

usage() {
	cat <<USAGE
Usage: $(basename "$0") check            # 🧪 fetch ${TAG} files into a temp repo, git apply --check every patch
       $(basename "$0") apply <CAMROOT>  # 🩹 git apply the patch series in a CAM7 checkout (idempotent)
USAGE
}

check() {
	# 🧪🌐 Sparse copy of the touched files at ${TAG}, then a dry run of the whole series.
	local f
	TMPD="$(mktemp -d "${TMPDIR:-/tmp}/cesm3_check.XXXXXX")"
	trap 'rm -rf "${TMPD}"' EXIT
	git -C "${TMPD}" init -q
	for f in "${FILES[@]}"; do
		mkdir -p "${TMPD}/$(dirname "${f}")"
		curl -sfL --max-time 60 "${RAW}/${f}" -o "${TMPD}/${f}" || { echo "ERROR: cannot fetch ${f}" >&2; exit 1; }
	done
	git -C "${TMPD}" apply --check -v "${HERE}"/patches/*.patch
	echo "OK: all patches apply to ${TAG}"
}

apply() {
	# 🩹🔁 Skip patches that are already in (reverse check succeeds), apply the rest.
	local cam="$1" p
	for p in "${HERE}"/patches/*.patch; do
		if git -C "${cam}" apply --reverse --check "${p}" 2>/dev/null; then
			echo "  already applied: $(basename "${p}")"
		else
			git -C "${cam}" apply "${p}"
			echo "  applied: $(basename "${p}")"
		fi
	done
}

case "${1:-}" in
check) check ;;
apply) [[ $# -eq 2 ]] || { usage; exit 1; }; apply "$2" ;;
*) usage; exit 1 ;;
esac
