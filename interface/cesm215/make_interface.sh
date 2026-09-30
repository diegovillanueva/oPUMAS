#!/usr/bin/env bash
# 🧊🔌 Glue for running oPUMAS inside CESM 2.1.5 CAM (tag cam_cesm2_1_rel_60).
# 📘📂 See README.md in this directory for the full workflow.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="${W:-/cluster/work/climate/dvillanueva/work/CESM}"   # CESM work dir
BASE="${BASE:-${W}/cesm215}"                              # existing CESM 2.1.5 root (source of externals)
TAG="${TAG:-cam_cesm2_1_rel_60}"                          # CAM tag the interface was written against
CAMWT="${CAMWT:-${W}/cam_rel60_pumas}"                    # CAM worktree
BRANCH="${BRANCH:-dv_pumas_rel60}"                        # branch of the CAM worktree
ROOT="${ROOT:-${W}/cesm215_pumas}"                        # sibling CESM root seen by cime

# 🧮📦 PUMAS files compiled by CAM (never micro_pumas_ccpp, micro_pumas_tempfix, micro_mg1_0).
PUMAS_FILES=(pumas_kinds pumas_gamma_function micro_pumas_utils micro_pumas_diags
	pumas_stochastic_collect_tau module_neural_net tau_neural_net_quantile
	ML_fixer_check micro_pumas_v1)
CAM_EXTERNALS=(chem_proc src/physics/clubb src/physics/carma/base src/physics/cosp2/src)
OVERRIDES=(micro_mg_cam.F90 wv_sat_methods.F90)
XML=bld/namelist_files/namelist_definition.xml

links() {
	# 🔗🧊 Relative symlinks to the oPUMAS sources (the only PUMAS files CAM sees).
	local f
	for f in "${PUMAS_FILES[@]}"; do
		ln -sfn "../../${f}.F90" "${HERE}/${f}.F90"
		[[ -r "${HERE}/${f}.F90" ]] || { echo "ERROR: missing ${f}.F90" >&2; exit 1; }
	done
	echo "${HERE}"
}

setup_tree() {
	# 🌿🧹 Clean CAM worktree at the stable tag, plus a sibling CESM root that points to it.
	local d c
	if [[ ! -d "${CAMWT}" ]]; then
		git -C "${BASE}/components/cam" worktree add -b "${BRANCH}" "${CAMWT}" "${TAG}"
	fi
	for d in "${CAM_EXTERNALS[@]}"; do
		[[ -e "${CAMWT}/${d}" ]] || ln -s "${BASE}/components/cam/${d}" "${CAMWT}/${d}"
	done
	mkdir -p "${ROOT}/components"
	for d in cime cime_config Externals.cfg describe_version; do
		[[ -e "${ROOT}/${d}" ]] || ln -s "${BASE}/${d}" "${ROOT}/${d}"
	done
	for c in "${BASE}"/components/*; do
		[[ "$(basename "${c}")" == cam ]] && continue
		[[ -e "${ROOT}/components/$(basename "${c}")" ]] || ln -s "${c}" "${ROOT}/components/"
	done
	[[ -e "${ROOT}/components/cam" ]] || ln -s "${CAMWT}" "${ROOT}/components/cam"
	echo "${CAMWT}"
	echo "${ROOT}"
}

patch_xml() {
	# 📝🧊 Add namelist group micro_pumas_nl to the worktree (idempotent), committed on BRANCH.
	local p="${HERE}/namelist_definition.micro_pumas_nl.patch"
	if grep -q 'group="micro_pumas_nl"' "${CAMWT}/${XML}"; then
		echo "already patched: ${CAMWT}/${XML}"
		return 0
	fi
	git -C "${CAMWT}" apply "${p}"
	git -C "${CAMWT}" commit -q -m "✨ feat: 🧊📝 add micro_pumas_nl namelist group for oPUMAS" -- "${XML}"
	echo "${CAMWT}/${XML}"
}

show_diff() {
	# 🔍📐 Effective patch of the two overridden CAM files against the tag.
	local f
	for f in "${OVERRIDES[@]}"; do
		diff -u --label "${TAG}:src/physics/cam/${f}" --label "interface/cesm215/${f}" \
			<(git -C "${CAMWT}" show "${TAG}:src/physics/cam/${f}") "${HERE}/${f}" || true
	done
}

case "${1:-}" in
	links) links ;;
	setup-tree) setup_tree ;;
	patch-xml) patch_xml ;;
	diff) show_diff ;;
	all) links; setup_tree; patch_xml ;;
	*)
		echo "Usage: $0 {all|links|setup-tree|patch-xml|diff}" >&2
		exit 1
		;;
esac
