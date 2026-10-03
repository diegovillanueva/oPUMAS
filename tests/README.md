# 🧪🧊 oPUMAS tests

🎯♻️ Goal: any commit can be checked for "still compiles" in one minute and for "still bit-for-bit" with one command.
🖥️⛔ Run everything inside `herdr-compute`, never on a login node.

## 🚀📋 Commands

```bash
tests/compile_check.sh            # 20 files (PUMAS core + CAM overrides), about 45 s, no Slurm job
tests/bfb.sh run                  # rebuild the reference case, run 3 time steps, compare cam.r with the baseline
tests/bfb.sh run sb2001           # the same with one switch from switches.csv set
tests/bfb.sh list                 # switch matrix and which baselines exist
tests/bfb.sh par sppe_            # all rows starting with sppe_ at once, one --keepexe clone per row (after 'run base')
STOP_OPTION=ndays STOP_N=3 WALLTIME=01:00:00 tests/bfb.sh par all   # 3-day check of every row (after 'run base' with the same settings)
python3 tests/cmp_restart.py a.cam.r.nc b.cam.r.nc   # compare any two restarts, also across cases
```

- 🧮✅ `compile_check.sh` reuses the compile line and `.mod` files of a built case (`CASE=`, default `cesm215_pumas.F2000climo.bfb_dev`). `IFACE=<worktree>/interface/cesm215` checks another branch. Right after `case.build --clean` the `.mod` files are gone, so it fails (`Cannot open module file`) until the case is built again.
- ♻️🧊 `bfb.sh` uses the case `cesm215_pumas.F2000climo.bfb_dev` (CESM 3-like namelist, graupel on, `-usr_src` = this repo). It submits one build job and one run job, so mind the 16-job limit.
- 📦📂 Baselines are restart files and live outside the repo: `/cluster/work/climate/$USER/cesm_runs/opumas_bfb_baselines/<case>/<STOP_OPTION><STOP_N>/<label>.cam.r.nc` (default `nsteps3`; 1-day baselines in `ndays1/`; 3-day baselines of every row, made before the single-`main` merge of 2026-10-03, in `ndays3/`).
- ⏱️🎚️ Run length: `STOP_OPTION=nsteps STOP_N=3` by default (1.5 model hours); set e.g. `STOP_OPTION=ndays STOP_N=1` for a longer check. The run job asks for `WALLTIME=00:30:00` so it backfills fast.
- 🐑⚡ `par <prefix>` (`par all` = every row except `base`) clones the reference case once per label into `/cluster/work/climate/$USER/cesm_runs/opumas_bfb_clones/<case>.<label>` (`--keepexe`: same executable and SourceMods, own run directory), runs `PAR=15` labels at a time and prints one result line per row; full logs in each clone's `bfb.log`. Rebuild the reference case before `par` after a Fortran change.
- 🔁📌 After an intended change of results, renew a baseline with `tests/bfb.sh save <label>`.
- ⚡🎚️ `BUILD=0 tests/bfb.sh run <label>` skips the build when only the namelist changed.
- 🧩🔁 Before each build `bfb.sh` deletes `bld/atm/obj/{Srcfiles,Depends}`: CIME never regenerates them, so a new source file (e.g. `micro_pumas_sppe.F90`) would otherwise fail with `Can't open module file`.

## 🎚️🧭 Strategy for testing switches

🪜📐 Four levels, from cheap to expensive; a commit climbs only as far as it needs.

| Level | Question | Tool | Cost | When |
|---|---|---|---|---|
| 0 | Does it compile? | `compile_check.sh` | 40 s | every commit |
| 1 | Is the default path unchanged? | `bfb.sh run` (label `base`) | build + a few min, 1 run job | every commit that touches Fortran |
| 2 | Is each switch wired, and unchanged? | `bfb.sh par <prefix>` over `switches.csv` | a few min plus queue, up to 15 run jobs | new switch; before and after a core refactor |
| 3 | Is the climate sensible? | 3-month runs in `interface/cesm215/TODO.md` | hours | new science feature |

- 🎛️🔕 **Level 1 is the golden rule as a command.** A commit that adds a switch must pass `base` without renewing the baseline, because the switch is off by default.
- 🔌✅ **Level 2 gives two checks per switch.** The restart must differ from `base`, which proves the namelist variable reaches the physics. It must also be identical to the switch's own baseline, which proves that later commits did not change that code path.
- 🧱🔎 **Why level 2 matters for refactoring.** The `base` namelist only exercises the branches it switches on. A refactor of `micro_pumas_tend` is safe only if every row of `switches.csv` stays bit-for-bit.
- ➕📄 **A new switch adds one row** to `switches.csv` (`label,expect,namelist`; several settings are separated by `;`, and `@REPO@` expands to the repo root). Settings are appended to the case's `user_nl_cam` for one run and removed afterwards.
- 🚦🔢 **Throttle.** `run` uses one job; `par` keeps at most `PAR=15` run jobs (16 with `herdr-compute`). Lower `PAR` when other jobs of yours are queued.

## 🧪⚡ Column driver

- 📦🔬 [`column/`](column/README.md) holds a standalone driver that runs `micro_pumas_tend` on 200 stored CAM columns in seconds. It reproduces CAM to 2.4e-16, writes one golden file per switch row, and records every process call for the Julia twin.

## 🚧📝 Not covered yet

- 🧊🚫 Switches that change the tracer count (graupel off) need a second reference case built without `-nadv 35`.
- 🤖🚫 `emulated` warm rain is not wired in the interface.
- 🪦🔕 `evap_sed_off` and `icenuc_use_meyers` have `expect=same`: upstream reads the namelist flags (`micro_pumas_v1.F90`, `micro_pumas_init`) but never uses them, so the 3-day run is identical to `base`.
