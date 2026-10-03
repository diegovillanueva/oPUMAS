# 🐙🔁 oPUMAS git workflow

🧭📌 `upstream` = ESCOMP/PUMAS (fetch only), `origin` = diegovillanueva/oPUMAS (private).
🌿🪞 `main_cam` mirrors `upstream/main_cam` exactly; all work happens on `main` (or `feat/*` branches off `main`).

## 🤖⏰ Agent routine: check upstream at every session start

🔍📥 Before any work in this repo, run:

```
git fetch upstream --tags
git log --oneline main..upstream/main_cam
```

- ✅🟢 Empty output: `main` already contains all upstream fixes, continue.
- 🔔🟡 Non-empty: tell Diego which upstream commits are new, then run the daily sync below.
- 🛑⚠️ If the merge into `main` conflicts, stop, show the conflicting files, and ask Diego before resolving.
- 🟣🔁 After a merge, run `julia/tools/check_twins.sh <old>..<new>`. Each printed line names a Fortran routine that changed and its oPUMAS.jl twin to update and re-test.
- 🔤🟣 After a merge or a dictionary edit, run `julia/tools/julia.sh --project=julia julia/tools/check_names.jl` (on herdr-compute): every `Julia` path must exist in oPUMAS.jl and every `Fortran` name in the sources.

## 🔄📅 Daily sync

```
git fetch upstream --tags
git switch main_cam && git merge --ff-only upstream/main_cam && git push origin main_cam --tags
git switch main && git merge main_cam
git push origin main
```

🔀🚫 Use merge, not rebase, on `main`: it is pushed, and rebasing would need force pushes (forbidden).

## ✍️🧩 Commit hygiene for future upstream PR

- 🧱✂️ Keep changes meant for ESCOMP in small, self-contained commits on `main`.
- 🙅📝 Do not mix oPUMAS-only files (`README.md`, `git_work.md`, `AGENTS.md`, handoffs, diagnostics) into those commits; they must never reach ESCOMP.

## 🔒⛔ Upstream PR: DISABLED

🚫🙋 **Do not create a public fork, add a `fork` remote, push to any ESCOMP-facing branch, or run `gh pr create` until Diego explicitly says so.**
📋🔮 Recipe for when he does (reference only):

```
gh repo fork ESCOMP/PUMAS --remote --remote-name fork      # one-time public fork
git fetch upstream
git switch -c pr/<topic> upstream/main_cam
git cherry-pick <commits from main>  # or squash into one clean commit
git push fork pr/<topic>
gh pr create --repo ESCOMP/PUMAS --head diegovillanueva:pr/<topic> --base main_cam
```

🔐🍴 Why a separate fork: GitHub forks of public repos cannot be private, and PRs must come from a fork, so `oPUMAS` (private, not a fork) cannot open the PR itself.
