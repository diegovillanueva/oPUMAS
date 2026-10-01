#!/usr/bin/env bash
# 🔁🌿 Two-way sync between the reference branch (first row of languages.csv) and every other language branch.
# 🧭📘 Usage: translation/sync.sh            (run from any worktree of the repository)
# 🧩🔗 The last synced pair is the newest commit on the language branch carrying a "Translated-from: <hash>" trailer.
# 1️⃣📥 Commits made on the language branch since then are replayed, translated, on top of the synced reference commit.
# 2️⃣🔀 They are merged into the reference branch (fast-forward when it has no new commits).
# 3️⃣📤 The reference tip is translated and committed on the language branch, so translate(reference) == language branch.
# 🚫📡 Nothing is pushed; push both branches to origin afterwards.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
top="$(git -C "${here}" rev-parse --show-toplevel)"
git_common="$(cd "${top}" && cd "$(git rev-parse --git-common-dir)" && pwd)"
tmp="${git_common}/translation-sync-tmp"

cleanup() { git -C "${top}" worktree remove --force "${tmp}" >/dev/null 2>&1 || true; }
trap cleanup EXIT

worktree_of() {  # path of the worktree that has branch $1 checked out, or empty
    git -C "${top}" worktree list --porcelain | awk -v b="refs/heads/$1" '/^worktree /{w=substr($0,10)} $1=="branch" && $2==b {print w}'
}

translated_tree() {  # $1 commit, $2 from-language, $3 to-language -> tree hash of the translated commit
    cleanup
    git -C "${top}" worktree add --quiet --detach "${tmp}" "$1"
    python3 "${tmp}/translation/translate.py" --from "$2" --to "$3" --root "${tmp}" >/dev/null
    git -C "${tmp}" add -A
    git -C "${tmp}" write-tree
    cleanup
}

advance() {  # move branch $1 to commit $2 (must be a fast-forward), also in the worktree that has it checked out
    local wt; wt="$(worktree_of "$1")"
    if [[ -n "${wt}" ]]; then
        git -C "${wt}" merge --quiet --ff-only "$2"
    else
        git -C "${top}" fetch --quiet . "$2:refs/heads/$1"
    fi
}

require_clean() {
    local wt; wt="$(worktree_of "$1")"
    if [[ -n "${wt}" && -n "$(git -C "${wt}" status --porcelain --untracked-files=no)" ]]; then
        echo "ERROR: branch $1 has uncommitted changes in ${wt}; commit them first." >&2; exit 1
    fi
}

mapfile -t langs < <(tail -n +2 "${here}/languages.csv" | cut -d, -f1,2)
ref_lang="${langs[0]%%,*}"; ref="${langs[0]##*,}"
require_clean "${ref}"

for entry in "${langs[@]:1}"; do
    lang="${entry%%,*}"; br="${entry##*,}"
    ref_tip="$(git -C "${top}" rev-parse "${ref}")"

    if ! git -C "${top}" rev-parse --verify --quiet "refs/heads/${br}" >/dev/null; then
        tree="$(translated_tree "${ref_tip}" "${ref_lang}" "${lang}")"
        new="$(printf '🔁 sync: create %s from %s %s\n\nTranslated-from: %s\n' "${br}" "${ref}" "${ref_tip:0:7}" "${ref_tip}" \
               | git -C "${top}" commit-tree "${tree}" -p "${ref_tip}")"
        git -C "${top}" branch "${br}" "${new}"
        echo "created ${br} at ${new:0:7} from ${ref} ${ref_tip:0:7}"
        continue
    fi

    require_clean "${br}"
    marker="$(git -C "${top}" log --first-parent -1 --format=%H --grep='^Translated-from: ' "${br}")"
    synced="$(git -C "${top}" log -1 --format=%B "${marker}" | sed -n 's/^Translated-from: //p' | tail -1)"
    br_tip="$(git -C "${top}" rev-parse "${br}")"
    if [[ "${marker}" == "${br_tip}" && "${synced}" == "${ref_tip}" ]]; then
        echo "${ref} and ${br} are already in sync"
        continue
    fi

    # 1. replay the new commits of the language branch, translated, on top of the synced reference commit
    tip="${synced}"
    for c in $(git -C "${top}" rev-list --reverse --first-parent "${marker}..${br}"); do
        tree="$(translated_tree "${c}" "${lang}" "${ref_lang}")"
        tip="$(git -C "${top}" log -1 --format=%B "${c}" | { cat; printf '\nTranslated-from: %s\n' "${c}"; } \
               | GIT_AUTHOR_NAME="$(git -C "${top}" log -1 --format=%an "${c}")" \
                 GIT_AUTHOR_EMAIL="$(git -C "${top}" log -1 --format=%ae "${c}")" \
                 GIT_AUTHOR_DATE="$(git -C "${top}" log -1 --format=%aI "${c}")" \
                 git -C "${top}" commit-tree "${tree}" -p "${tip}")"
        echo "replayed ${br} ${c:0:7} as ${tip:0:7}"
    done

    # 2. bring them into the reference branch
    if [[ "${tip}" != "${synced}" ]]; then
        if git -C "${top}" merge-base --is-ancestor "${ref_tip}" "${tip}"; then
            advance "${ref}" "${tip}"
        else
            wt="$(worktree_of "${ref}")"
            [[ -n "${wt}" ]] || { echo "ERROR: check out ${ref} in a worktree to merge." >&2; exit 1; }
            git -C "${wt}" merge --no-edit -m "🔀 merge: commits from ${br} (translated)" "${tip}" \
                || { echo "ERROR: resolve the merge in ${wt}, commit, and run sync.sh again." >&2; exit 1; }
        fi
        ref_tip="$(git -C "${top}" rev-parse "${ref}")"
    fi

    # 3. translate the reference tip onto the language branch (this commit is the new sync marker)
    tree="$(translated_tree "${ref_tip}" "${ref_lang}" "${lang}")"
    new="$(printf '🔁 sync: %s %s\n\n%s\n\nTranslated-from: %s\n' "${ref}" "${ref_tip:0:7}" \
           "$(git -C "${top}" log --format='- %h %s' "${synced}..${ref_tip}")" "${ref_tip}" \
           | git -C "${top}" commit-tree "${tree}" -p "${br_tip}")"
    advance "${br}" "${new}"
    echo "synced ${ref} ${ref_tip:0:7} -> ${br} ${new:0:7}"
done
echo "next: git push origin ${ref} $(printf '%s ' "${langs[@]:1}" | sed 's/[^ ]*,//g')"
