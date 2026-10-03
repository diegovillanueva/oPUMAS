# 🤝🧊 Handoff: Sppe controls on oPUMAS (`dev/Sppe`), switch rows verified

## 📌 Commit
`bc15c10` ✨ feat: ⏱️🐑 bfb.sh runs 3 time steps by default and switch rows in parallel (par)
(🔗 previous: `handoffs/sppe.949ac2e.md`, the Sppe implementation itself)

## ✅ What was done
- ⏱️🧪 `tests/bfb.sh` runs 3 time steps by default (`STOP_OPTION`/`STOP_N`), with baselines per length in `opumas_bfb_baselines/<case>/nsteps3/`; the 1-day ones are in `ndays1/`.
- 🐑⚡ `tests/bfb.sh par sppe_` runs all `sppe_*` rows at once, one `--keepexe` clone each in `/cluster/work/climate/dvillanueva/cesm_runs/opumas_bfb_clones`, at most 15 run jobs.

## 🧪 Current state
- ✅ `base` (nsteps3) stored from the `dev/Sppe` build; the 1-day `base` had already passed bit-for-bit, so the defaults are neutral.
- ✅ 15 `sppe_*` rows differ from `base` (each switch reaches the physics); `sppe_no_cnt_frz` is identical to `base`, as expected with `use_hetfrz_classnuc = .true.`. All 16 baselines stored; about 15 min wall time.
- ✅ `tests/compile_check.sh` OK.

## ➡️ Next steps
- 🧪 One run with `interface/cesm215/user_nl_cam.pumas_sppe`; check the `micro_dv_nl` values in `atm.log`.
- 🐙 Sppe_revision_public: commit the local `PUMAS` branch (scratchpad clone), rename `main` to `cesm215_mg2`, push `PUMAS` and make it the default.
- 🧊 Open in `interface/cesm215/TODO.md`: graupel riming/HM, `rainfreeze_ifs`, TAU.
- 🌍 Check whether CAM7's `micro_pumas_cam.F90` compiles against oPUMAS.

## 🧭 Key decisions
- ⏱️ "3 day test" read as 3 time steps (the test was already 1 day and the request was "faster"); `STOP_OPTION=ndays STOP_N=3` gives 3 days.
- 📂 `create_clone` resets `SRCROOT` to the stock tree; `par` rewrites it to the reference case's CAM worktree.
- 🟰 `same` rows are now checked against `base` too.
