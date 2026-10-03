# 🌿1️⃣ Handoff: one development branch `main`

## 📌 Commit
- 🔀 `4178f90` merge of `dev/Sppe`, `1766338` column driver always builds Sppe (3-day loop), `1685607` translator retired; this handoff is the commit after them.

## ✅ What was done
- 🌿🏷️ `main_cam_interface` was renamed to `main` (local, `origin`, GitHub default branch). `dev/Sppe` and `dev/julia` were merged into it.
- 🗄️🏷️ Old branches are now the tags `archive/dev-julia`, `archive/dev-Sppe` and `archive/CF-CCPP`. The branches were deleted locally and on `origin`, and the worktree `cesm215_pumas_worktrees/opumas_CF-CCPP` was removed.
- 🔤📖 The CF-CCPP translator is gone. `translation/dictionary.csv` (352 rows, columns `id,kind,scope,units,meaning,Fortran,Julia,standard_name,notes`) maps upstream Fortran names to the Julia twin, and test set M7 of `Pkg.test()` (`julia/test/names_check.jl`) checks both directions.

## 🧪 Validation (3-day runs, `STOP_OPTION=ndays STOP_N=3`)
- 🧊 Before the merge: baselines in `cesm_runs/opumas_bfb_baselines/cesm215_pumas.F2000climo.bfb_dev/ndays3/`, from `main_cam_interface` (base and the CESM rows) and from `dev/Sppe` (base bit-for-bit with the first, plus the `sppe_*` rows).
- ✅ After the merge: base and 31 rows PASS bit-for-bit against those baselines. `evap_sed_off` and `icenuc_use_meyers` are identical to base (`expect=same`: these upstream flags are never used).
- ✅ `tests/compile_check.sh` OK (20 files), and `Pkg.test()` 7699/7699, including the column loop of 432 steps (3 days at 600 s), bit-for-bit 155/155.

## ➡️ Next steps
- 🧊 Open items stay in `interface/cesm215/TODO.md` and `handoffs/julia.plan.md`.
- 🗑️ The CESM case `cesm215_pumas.F2000climo.bfb_cfccpp` pointed at the removed CF-CCPP worktree; delete it or restore it from `archive/CF-CCPP` if it is still needed.

## 🧭 Key decisions
- 🔤 Readable names live only in the Julia twin, and the Fortran keeps the upstream names, so upstream merges stay simple.
- 🧪 The 3-day baselines (`ndays3/`) sit next to the 3-step ones (`nsteps3/`); `bfb.sh` chooses the folder from `STOP_OPTION`/`STOP_N`.
- 🧩 `tests/compile_check.sh` reuses the `.mod` files of the built case, so it fails right after `case.build --clean` until the case is rebuilt.
