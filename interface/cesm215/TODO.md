# 📝🧭 Open work for the CESM 2.1.5 interface

📌🔄 Single source of truth for what the interface does not do yet. Update it in the commit that implements or tests an item; handoffs (`../../handoffs/`) only point here.
🎛️🔕 Every item must follow the rules in `../../AGENTS.md`: switchable, off by default, CESM 2.3 defaults.

## 🧪🎚️ Deferred: fill the switch baselines

- 📌⏳ Status: open, to do later. Only the label `base` has a baseline (2026-10-01).
- 🎯 Goal: every row of `../../tests/switches.csv` has run once, so each switch is proven to be wired and its code path is protected bit-for-bit.
- ▶️ Proposed test, one label at a time (about 5 minutes each, one run job, no rebuild):
  ```bash
  cd /cluster/work/climate/dvillanueva/work/CESM/PUMAS
  for l in $(tail -n +3 tests/switches.csv | cut -d, -f1); do BUILD=0 tests/bfb.sh run "$l" || echo "FAILED $l"; done
  tests/bfb.sh list
  ```
- ✅ Pass criterion per label: the run ends, the restart differs from `base`, and the baseline is stored (`NEW`). A second run of the same label must report `PASS`.
- 🔧 Expected follow-up: fix rows whose namelist line is rejected or has no effect, then move the label from "Wired but not yet run" below to tested.
- 🚧 Not covered: graupel off needs a second reference case without `-nadv 35`.

## 💡🧩 Deferred: replace the full-file CAM overrides

- 📌⏳ Status: open, to do later. Full plan and trade-offs: `PROPOSAL_overrides.md`.
- 🪜 Plan:
  1. Fill the switch baselines (item above), so every code path is protected.
  2. Store the CAM edits as a patch series in `patches/` and apply them with `git am` to the worktree branch `dv_pumas_rel60`; add `export` and `verify` to `make_interface.sh`; remove the five full copies. Test: `tests/bfb.sh run` passes without a new baseline.
  3. Move the PUMAS-specific code of `micro_mg_cam.F90` into a new module `micro_pumas_cam.F90` in small moves (readnl, init, tend). Test after each move: `tests/compile_check.sh` and `tests/bfb.sh run`.
- ⚠️ Risk to settle in step 3: the glue uses module-private variables of `micro_mg_cam` (pbuf and constituent indices).

## 🔧🧊 Features that need interface changes

| Feature | What is missing | Effort | Status |
|---|---|---|---|
| Emulated (neural-net) warm rain, `micro_pumas_warm_rain='emulated'` | Three emulator file names are passed as `''` to `micro_pumas_init`; readnl rejects `emulated`; the NetCDF emulator files are not on disk | small code, files to obtain | open |
| Constant rain/snow number (`nrcons`, `nrnst`, `nscons`, `nsnst`) | Hard-wired `.false., 0._r8` in the init call; needs four namelist variables and XML entries | small | open |
| Radiatively active graupel (`graupel_in_rad`) | pbuf `DEGRAU`, `ICGRAUWP`, `CLDFGRAU` in `micro_mg_cam.F90`; overrides of rel_60 `cloud_rad_props.F90` and `radiation.F90` (port from cam6_3_097). Off by default upstream too | medium | open |
| Graupel/hail diagnostics | `LS_REFFGRAU`, `CLDFGRAU`, `UMG`, `FREQG`, `QGSEDTEN`, graupel process rates (`PRACGO`, `PSACWGO`, `PGSACWO`, `PGRACSO`, `PRDGO`, `QMULTGO`, `QMULTRGO`, `PSACRO`, `NPRACGO`, `NSCNGO`, `NGRACSO`, `NMULTGO`, `NMULTRGO`, `NPSACWGO`, `MELTGTOT`); now local dummies | small to medium | open |
| PUMAS-only process rates and TAU diagnostics | `vapdepstot`, `mnudeptot`, `mnuccritot`, `nraggtot`, number rates (`n*tot`) and the bin-scheme tendencies exist in `proc_rates` but are not written | small to medium | open |
| 10 cm radar reflectivity (`refl10cm`, `reflz10cm`) | Computed by PUMAS, discarded | small | open |
| COSP sees graupel | rel_60 `cospsimulator_intr` knows only rain and snow; the graupel flux is folded into `LS_FLXSNW` | medium | open |

## 🧊🎛️ Sppe controls on PUMAS-only paths (`micro_pumas_sppe.F90`)

- 🧊❄️ Graupel riming (`psacwg`, `pgsacw`) is not scaled by `riming_factor`, and graupel Hallett–Mossop (`qmultg`, `qmultrg`) ignores `onoff_sip`/`factor_HM`; the MG2 patches never had graupel.
- 🌧️🧊 `micro_mg_rainfreeze_ifs` freezes rain outside `heterogeneous_rain_freezing`, so `rain_freeze = .false.` does not stop it.
- 🔢🌧️ `auto_liq_factor` does not touch the TAU bin warm rain.

## 🧪⏳ Wired but not yet run in CESM 2.1.5

| Feature | Namelist | Status |
|---|---|---|
| Hail instead of graupel | `micro_mg_do_hail` | untested |
| Constant graupel/hail number | `micro_mg_ngcons`, `micro_mg_ngnst` | untested |
| Seifert–Beheng warm rain | `micro_pumas_warm_rain='sb2001'` | untested |
| TAU bin warm rain | `micro_pumas_warm_rain='tau'`, `micro_pumas_kernel_file` (`KBARF_tau_kernel.dat` in the repo root) | untested |
| IFS-style options | `micro_mg_evap_scl_ifs`, `micro_mg_evap_rhthrsh_ifs`, `micro_mg_rainfreeze_ifs`, `micro_mg_ifs_sed` | untested |
| Evaporation and ice nucleation variants | `micro_mg_evap_sed_off`, `micro_mg_icenuc_rh_off`, `micro_mg_icenuc_use_meyers` | untested |
| Constant droplet/ice number | `micro_mg_nccons`, `micro_mg_nicons` | untested |
| CESM 3-like set | `user_nl_cam.pumas_cesm3like` (implicit fall, `accre_sees_auto`, `precip_fall_corr`, CAM7 tuning, `microp_aero_wsub_*`) | 1 day OK (case `pumas_c3like`, 2026-10-01); multi-month run open |

## 🌍📊 Runs and comparisons

- 3-month CESM 3-like run (`pumas_c3like`) and comparison with `pumas_test`.
- 3-month graupel run (`pumas_grau`) and comparison with `pumas_test` (graupel off).
- `pumas_test` against an MG2 reference (`PRECT`, `CLDLIQ`, `CLDICE`, `TGCLDLWP`, `SWCF`, `LWCF`).
