📘🔗 Follow `/cluster/home/dvillanueva/config_and_bash/AGENT_INSTRUCTIONS.md`.
🔄🤖 Before any work, run the upstream check in `git_work.md` ("Agent routine") and report new ESCOMP/PUMAS commits.
⛔🐙 Never push to `upstream` (ESCOMP/PUMAS); push only to `origin` (diegovillanueva/oPUMAS).
🔒🚫 Upstream PRs are disabled: no public fork, no `fork` remote, no `gh pr create` until Diego explicitly says so.
🪞🌿 Never commit on `main_cam`; it is a fast-forward mirror of `upstream/main_cam`.
🎛️🔕 Golden rule: every new feature built on oPUMAS is switchable and off by default, so an empty namelist makes oPUMAS behave exactly like upstream PUMAS (bit-for-bit where possible). Add a namelist flag (default `.false.`/neutral value), wire it through the host interface, and document it in `interface/*/README.md`.
🧭🎚️ Default-behaviour rule: every added feature, in oPUMAS core or in any `interface/*`, is switchable and by default behaves like upstream PUMAS. When out-of-the-box CESM–PUMAS versions differ (e.g. `micro_mg_do_graupel=.true.` in CESM 3 / CAM7 but `.false.` in CESM 2.3 / cam6), keep the **CESM 2.3 behaviour as the default for now**; the CESM 2.1.5 interface follows the CESM 2.3 defaults.
