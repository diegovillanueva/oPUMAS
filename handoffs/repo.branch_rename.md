# 🌿🏷️ Handoff: `dev` renamed to `main_cam_interface`

- 📌 Commit: see `git log -1 -- handoffs/repo.branch_rename.md`.
- 🔁 What was done: the development branch `dev` was renamed to `main_cam_interface` (local, `origin`, GitHub default branch). The live docs and `translation/languages.csv` now point to the new name, and `sync.sh` carried the change to `CF-CCPP`.
- 📍 Current state: `main_cam_interface` = `main_cam` + interface, translation, tests and docs; the upstream PUMAS sources are still identical to `main_cam`.
- 🧭 Key decisions: older handoffs and the dated `code_analysis/*.md` notes keep the name `dev` as history. The CESM test case `cesm215_pumas.F2000climo.bfb_dev` keeps its name, so existing baselines stay valid.
