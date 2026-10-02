#!/usr/bin/env bash
# 🔗🔄 List the Julia twins whose Fortran routine changed in a commit range.
#    Usage: julia/tools/check_twins.sh <old>..<new>   (e.g. after an upstream merge)
#    Reads julia/TWINS.csv; prints "routine (file) -> julia_function" per change.
#    Variants such as _vect, _line, _2d and _r8 count as the routine itself.
set -euo pipefail

range="${1:?usage: check_twins.sh <old>..<new>}"
old="${range%%..*}"
new="${range##*..}"
new="${new:-HEAD}"
repo="$(git rev-parse --show-toplevel)"
twins="${repo}/julia/TWINS.csv"

# 📍 Path of a Fortran file named in TWINS.csv: repo root, else the CESM 2.1.5 interface.
fortran_path() {
    local name="$1" rev="$2"
    for path in "${name}" "interface/cesm215/${name}"; do
        if git -C "${repo}" cat-file -e "${rev}:${path}" 2>/dev/null; then
            echo "${path}"
            return
        fi
    done
}

# 🗺️ "start end name" of every subroutine and function of a file at a revision.
routine_ranges() {
    git -C "${repo}" show "$1:$2" | awk '
        tolower($0) ~ /^[ \t]*end[ \t]+(subroutine|function)/ {
            if (name != "") print start, NR, name; name = ""; next }
        BEGIN { prefix = "^[ \t]*((pure|elemental|recursive)[ \t]+)*"
                start_re = prefix "(subroutine|function)[ \t]" }
        tolower($0) ~ start_re {
            line = $0; sub(/\(.*/, "", line); n = split(line, w, /[ \t]+/)
            name = w[n]; start = NR }'
}

# 🔢 Changed line numbers on one side ("old" or "new") of the diff of a file.
changed_lines() {
    git -C "${repo}" diff -U0 "${old}" "${new}" -- "$1" | awk -v side="$2" '
        /^@@/ {
            spec = (side == "old") ? $2 : $3; sub(/^[-+]/, "", spec)
            n = split(spec, a, ","); start = a[1]; count = (n > 1) ? a[2] : 1
            if (count == 0) print start; else for (i = 0; i < count; i++) print start + i }'
}

# 🧮 Routines of a file touched by the changed lines on one side.
touched_routines() {
    local path="$1" rev="$2" side="$3"
    awk 'NR == FNR { s[NR] = $1; e[NR] = $2; r[NR] = $3; n = NR; next }
         { for (i = 1; i <= n; i++) if ($1 >= s[i] && $1 <= e[i]) print r[i] }' \
        <(routine_ranges "${rev}" "${path}") <(changed_lines "${path}" "${side}")
}

found=0
for file in $(tail -n +2 "${twins}" | cut -d, -f2 | sort -u); do
    [[ "${file}" == *.F90 ]] || continue
    path_new="$(fortran_path "${file}" "${new}")"
    path_old="$(fortran_path "${file}" "${old}")"
    routines="$( (
        if [[ -n "${path_new}" ]]; then touched_routines "${path_new}" "${new}" new; fi
        if [[ -n "${path_old}" ]]; then touched_routines "${path_old}" "${old}" old; fi
    ) | tr 'A-Z' 'a-z' | sort -u)"
    for routine in ${routines}; do
        twin_rows="$(awk -F, -v f="${file}" -v r="${routine}" \
            'BEGIN { base = r
                     variants = "line|vect2|vect|vec|2d|r8_vec|r8_vect|r8"
                     variants = variants "|integer_vect|integer"
                     sub("_(" variants ")$", "", base) }
             NR > 1 && $2 == f && (tolower($3) == r || tolower($3) == base) { print $1 }' \
            "${twins}" | sort -u)"
        for twin in ${twin_rows}; do
            echo "🔁 ${routine} (${file}) -> ${twin}"
            found=1
        done
    done
done
[[ "${found}" == 1 ]] || echo "✅ no translated routine changed in ${range}"
