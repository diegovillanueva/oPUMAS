📘🔗 Follow `/cluster/home/dvillanueva/config_and_bash/AGENT_INSTRUCTIONS.md`.
🔄🤖 Before any work, run the upstream check in `git_work.md` ("Agent routine") and report new ESCOMP/PUMAS commits.
⛔🐙 Never push to `upstream` (ESCOMP/PUMAS); push only to `origin` (diegovillanueva/oPUMAS).
🔒🚫 Upstream PRs are disabled: no public fork, no `fork` remote, no `gh pr create` until Diego explicitly says so.
🪞🌿 Never commit on `main_cam`; it is a fast-forward mirror of `upstream/main_cam`.
🎛️🔕 Golden rule: every new feature built on oPUMAS is switchable and off by default, so an empty namelist makes oPUMAS behave exactly like upstream PUMAS (bit-for-bit where possible). Add a namelist flag (default `.false.`/neutral value), wire it through the host interface, and document it in `interface/*/README.md`.
🧭🎚️ Default-behaviour rule: every added feature, in oPUMAS core or in any `interface/*`, is switchable and by default behaves like upstream PUMAS. When out-of-the-box CESM–PUMAS versions differ (e.g. `micro_mg_do_graupel=.true.` in CESM 3 / CAM7 but `.false.` in CESM 2.3 / cam6), keep the **CESM 2.3 behaviour as the default for now**; the CESM 2.1.5 interface follows the CESM 2.3 defaults.
🔤📖 Naming rule: follow `PHILOSOPHY.md`. Names are changed only through `translation/dictionary.csv` on `dev` (never by hand in the code), checked with `python3 translation/check.py`, and never touch namelist names, history field names, module names or the `micro_pumas_init`/`micro_pumas_tend` interfaces.
🌿🔁 Two editable branches: `dev` (PUMAS names) and `CF-CCPP` (translated). Run `translation/sync.sh` before and after working on `CF-CCPP`, then `python3 translation/check.py --branches`, then push both to `origin`. Upstream merges and dictionary edits go to `dev` only.
🤝📂 Handoffs live in `handoffs/<topic>.<shorthash or name>.md` at the repo root, not next to the scripts (this overrides the global rule).
🧪♻️ Before committing Fortran changes run `tests/compile_check.sh`; before committing anything that could change results run `tests/bfb.sh run` (see `tests/README.md`).
