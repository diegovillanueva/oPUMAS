# 🤝🧊 Handoff: Sppe controls on oPUMAS (`dev/Sppe`)

## 📌 Commit
`949ac2e` ✨ feat: 🧊🎛️ Sppe process controls on oPUMAS (micro_pumas_sppe) with CESM 2.1.5 wiring and CESM 3 patches

## ✅ What was done
- 🧊🎛️ The 20 Sppe controls (group `micro_dv_nl`, from github.com/diegovillanueva/Sppe_revision_public) work on PUMAS: 11 in the oPUMAS core via the new module `micro_pumas_sppe.F90` and setter `micro_pumas_sppe_set`, 9 CAM-side.
- 🔌🧊 CESM 2.1.5: `interface/cesm215` reads the namelist (`micro_dv_nml.F90`, hook in `microp_aero`), hands it to PUMAS from `micro_mg_cam_init`, and overrides `hetfrz_classnuc(_cam)`, `clubb_intr`, `macrop_driver`. Full set: `interface/cesm215/user_nl_cam.pumas_sppe`.
- 🌍🩹 CESM 3 / CAM7: `interface/cesm3/patches` against `cam6_4_204`; `make_interface.sh check` passes, `apply` is idempotent. Not compiled.

## 🧪 Current state
- ✅ `tests/bfb.sh run` (base): PASS, 236 variables bit-for-bit, so defaults are neutral.
- ✅ `tests/compile_check.sh` OK, `translation/check.py` 0 errors.
- ⏳ The `sppe_*` rows of `tests/switches.csv` were running when this was written (results in the follow-up commit).

## ➡️ Next steps
- 🧪 Finish the `sppe_*` switch rows; one 1-day run with `user_nl_cam.pumas_sppe` (check the `micro_dv_nl` values in `atm.log`).
- 🐙 Sppe_revision_public: rename `main` to `cesm215_mg2`, new default branch `PUMAS` (patch sets `pumas/`, `cam_cesm2_1_rel_60/`, `cam6_4_204/`).
- 🧊 Open in `interface/cesm215/TODO.md`: graupel riming/HM, `rainfreeze_ifs`, TAU are not reached by the controls.
- 🌍 Check whether CAM7's `micro_pumas_cam.F90` compiles against oPUMAS (CAM7 pins PUMAS v1.17).

## 🧭 Key decisions
- 🔌 New module + setter, not new `micro_pumas_init` arguments, so the PUMAS interfaces stay upstream-compatible.
- 🌿 `dev/Sppe` is not synced to `CF-CCPP`; namelist names stay as in Sppe (`micro_dv_nl`), so the Sppe `user_nl_cam` works unchanged.
- 🌍 In CAM7, `limfacdu` is read but unused; the native `hetfrz_dust_scalfac` replaces it.
- 🧩 `bfb.sh` deletes CIME `Srcfiles`/`Depends` before a build (stale lists broke the first build).
