# 🤝🔤 Handoff: name dictionary, translation scripts and the CF-CCPP branch (2026-10-01)

- 🔖 Commits on `dev`: `8ebc97e` (dictionary, scripts, docs), `6797ee5` (fix for components of other derived types), plus the docs commit that holds this file. `CF-CCPP` is created and updated only by `translation/sync.sh`.
- ✅ What was done: `translation/` with `dictionary.csv` (219 names), `translate.py`, `check.py`, `sync.sh`, tests and README; branch `CF-CCPP`; worktree `/cluster/work/climate/dvillanueva/work/CESM/cesm215_pumas_worktrees/opumas_CF-CCPP`; CESM cases `cesm215_pumas.F2000climo.bfb_dev` and `cesm215_pumas.F2000climo.bfb_cfccpp`.
- 📍 Current state: pytest (11 tests) and `check.py` pass; `check.py --branches` passes; `sync.sh` was exercised in a throwaway clone (dev-only commit, CF-CCPP-only commit, commits on both sides with a merge, nothing to do). CESM 2.1.5 bit-for-bit test passed: after 1 day (graupel plus CESM 3-like namelist) all 236 data variables of the `cam.r` restart of `bfb_cfccpp` are bit-identical to `bfb_dev`; only the path string `hrestpath` and the attribute `caseid` differ.
- ➡️ Next steps: Diego reviews the names in `dictionary.csv` (column `CF-CCPP`); next batches are scratch variables, TAU bin arrays, and the remaining about 600 identifiers; an `AMS` language is one new column plus one row in `languages.csv`; routine-level scope (`file.F90:routine`) is not implemented, only `global`, `type:` and `file:`.
- 🧠 Key decisions:
  - 🌿 `sync.sh` writes one squashed `🔁 sync` commit on `CF-CCPP` per sync instead of replaying each `dev` commit; commits made on `CF-CCPP` are replayed one by one onto `dev`.
  - 🏷️ Local in-cloud/in-precipitation rates carry the prefix `local_`; `proc_rates` components (grid means) carry the plain name.
  - 🔤 `cldliq` and `cldice` are abbreviated so that every name plus unit suffix fits in 63 characters (the longest name has exactly 63).
  - 🧱 Inside a derived-type definition only that type's own dictionary components are translated (`MGHydrometeorProps%rho` stays).
  - 🚫 `dv` is left out because the source also spells it `Dv`.
  - 🗑️ `handoff.philosophy.md` and `name_map.csv` are superseded by this directory.
