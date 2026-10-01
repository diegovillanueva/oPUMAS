# 🧊🔌 oPUMAS in CESM 2.1.5

🎯🧊 This directory lets CESM 2.1.5 CAM (tag `cam_cesm2_1_rel_60`) run the current PUMAS (`micro_pumas_v1`) as `micro_mg_version = 3`.
🧹📦 The oPUMAS sources are untouched, and CAM needs one namelist XML patch.

## 🚀📋 Workflow

```bash
W=/cluster/work/climate/dvillanueva/work/CESM
I=$W/PUMAS/interface/cesm215

# once: symlinks, clean CAM worktree at the tag, sibling CESM root, namelist XML patch
$I/make_interface.sh all

# every case: the usual wrapper, pointed at the sibling root and this dir
cd $W
cesmvvv=cesm215_pumas_worktrees/cesm215_pumas CAM_USR_SRC=$I \
  ./prepare_cesm.sh create my_case F2000climo --cam-namelist=$I/user_nl_cam.pumas
./prepare_cesm.sh submit cesm215_pumas.F2000climo.my_case   # with cesmvvv=cesm215_pumas_worktrees/cesm215_pumas too
```

🔍✅ After `create`, check that the case sees the right tree:

```bash
./xmlquery SRCROOT,CAM_CONFIG_OPTS   # SRCROOT = $W/cesm215_pumas_worktrees/cesm215_pumas, CAM_CONFIG_OPTS ends in -usr_src $I
```

## 🧩⚙️ How it works

- 🗂️🔝 CAM `configure -usr_src <dir>` puts `<dir>` first in `Buildconf/camconf/Filepath`, so a file here replaces the CAM file of the same name.
- 🔗🧊 The nine PUMAS files CAM compiles are relative symlinks to the repo root: `pumas_kinds`, `pumas_gamma_function`, `micro_pumas_utils`, `micro_pumas_diags`, `pumas_stochastic_collect_tau`, `module_neural_net`, `tau_neural_net_quantile`, `ML_fixer_check`, `micro_pumas_v1`.
- 🚫📄 `micro_pumas_ccpp`, `micro_pumas_tempfix` and the repo-root `micro_mg1_0` are never linked, because they need CCPP or clash with CAM.
- 🛠️📄 Five CAM files are overridden with full copies of the rel_60 file plus additions. `./make_interface.sh diff` prints the exact patch against the tag.
  - 🌧️🧊 `micro_mg_cam.F90` gains a `case (3)` in readnl, init, get_cols and tend_pack, and reads the extra namelist group `micro_pumas_nl`.
  - 💧📐 `wv_sat_methods.F90` gains `wv_sat_qsat_water_vect` and `wv_sat_qsat_ice_vect`, thin wrappers around the elemental rel_60 routines.
  - ⚖️❄️ `check_energy.F90` counts `GRAUQM` in the ice water budget and excludes it from the dry-mass fixer, so graupel conserves mass and energy.
  - 💨🌫️ `microp_aero.F90` and `ndrop.F90` expose the CAM7 sub-grid vertical velocity controls for droplet and ice activation (`microp_aero_wsub_*`). The defaults are the values rel_60 hard-codes, so they change nothing unless set.
- 📝🧊 `namelist_definition.micro_pumas_nl.patch` adds the `micro_pumas_nl` group and the `microp_aero_wsub_*` entries to CAM. It is committed on the worktree branch `dv_pumas_rel60`.
- 🌳🔗 `$W/cesm215_pumas_worktrees/cesm215_pumas` is a CESM root made of symlinks to `$W/cesm215`, except `components/cam`, which points to the clean worktree `$W/cesm215_pumas_worktrees/cam_rel60_pumas`.
- 🧭🐍 CIME resolves the `cime` symlink, so the wrapper passes `--srcroot` to `create_newcase`. Cases therefore sit physically in `$W/cesm215/cime/scripts` with the prefix `cesm215_pumas.`.

## 🎛️📝 Namelist

🔢🧊 `user_nl_cam.pumas` sets `micro_mg_version = 3`. Optional PUMAS knobs go in `user_nl_cam` as well; they default to the CAM6 values.

| variable | default | note |
|---|---|---|
| `micro_pumas_warm_rain` | `kk2000` | `kk2000`, `sb2001` or `tau` |
| `micro_pumas_kernel_file` | none | path to `KBARF_tau_kernel.dat`, needed for `tau` |
| `micro_mg_accre_enhan_fact`, `micro_mg_autocon_fact`, `micro_mg_autocon_nd_exp`, `micro_mg_autocon_lwp_exp` | 1, 0.01, -1.1, 2.47 | warm rain tuning |
| `micro_mg_homog_size`, `micro_mg_max_nicons` | 25e-6, 1e8 | ice |
| `micro_mg_vtrmi_factor`, `micro_mg_vtrms_factor`, `micro_mg_effi_factor`, `micro_mg_iaccr_factor` | 1 | ice and snow scaling |
| `micro_mg_evap_sed_off`, `micro_mg_icenuc_rh_off`, `micro_mg_icenuc_use_meyers`, `micro_mg_*_ifs`, `micro_mg_precip_fall_corr`, `micro_mg_implicit_fall`, `micro_mg_accre_sees_auto` | `.false.` | switches |
| `micro_mg_do_graupel`, `micro_mg_do_hail` | `.false.` | prognostic graupel or hail (one at most); adds constituents `GRAUQM`, `NUMGRA` |
| `micro_mg_ngcons`, `micro_mg_ngnst` | `.false.`, 0.0005e6 | constant graupel/hail number |
| `microp_aero_wsub_scale`, `microp_aero_wsubi_scale` | 1, 1 | scale factor on sub-grid vertical velocity for droplet activation and ice nucleation (group `microp_aero_nl`) |
| `microp_aero_wsub_min`, `microp_aero_wsub_min_asf`, `microp_aero_wsubi_min` | 0.2, 0.1, 0.001 m/s | minima before scaling, after scaling (in `dropmixnuc`), and for ice |

🧊❄️ `user_nl_cam.pumas_graupel` is the same plus `micro_mg_do_graupel = .true.` (10 constituents). 🔢⚠️ CAM's `configure` counts tracers before it sees this switch, so a graupel case also needs `./xmlchange --append CAM_CONFIG_OPTS="-nadv 35"` (the case's `nadv` in `Buildconf/camconf/config_cache.xml` plus 2) before `case.setup --reset` and a `case.build --clean-all` rebuild (make does not notice the changed `-DPCNST`, so an incremental build keeps the old executable); otherwise init aborts with `cnst_add: advected tracer index greater than pcnst`. 🌍🧊 `user_nl_cam.pumas_cesm3like` sets everything CESM 3 (CAM7) switches on out of the box: graupel, implicit fall, `accre_sees_auto`, `precip_fall_corr`, the CAM7 tuning factors and the `microp_aero_wsub_*` values. 📌🔖 Reference: ESCOMP/CESM `cesm3_0_beta10`, which pins ESCOMP/CAM `cam6_4_203` (microphysics defaults identical in `cam6_4_204`); CAM7 itself runs PUMAS `pumas_cam-release_v1.17_rename`, this repo is at v1.41+. It needs the same `-nadv 35` and clean rebuild as the graupel case. 🤖🚫 `emulated` warm rain is not wired. The existing MG2 namelist (`micro_mg_dcs`, `micro_mg_berg_eff_factor`, number constants, `micro_mg_num_steps`) is passed through unchanged.

## 🗺️🔀 Output mapping

- 📊✅ All standard MG2 history fields keep their names. PUMAS process rates in `proc_rates` are copied to the MG2 arrays, e.g. `pratot` to `PRAO` and `bergstot` to `BERGSO`.
- 🧊✂️ PUMAS splits snow deposition into `bergs` and `vapdeps`, so `BERGSO` is not comparable to MG2. PUMAS-only rates such as `vapdepstot`, `mnudeptot` and the number rates are not written.
- ❄️📊 With graupel on, `AQGRAU` and `ANGRAU` are written and the graupel flux is added to `LS_FLXPRC` and `LS_FLXSNW`. Graupel is not seen by radiation (no `ICGRAUWP`/`DEGRAU` fields, rel_60 `cloud_rad_props` does not know them).
- 🗑️📡 The 10 cm reflectivity, graupel effective radius and diameter go to local dummies.

## 🔎🧭 Audit against CESM 2.3 / CESM 3 (2026-10-01)

- 🧪📖 References: `cam6_3_097` (CESM 2.3-era, local copy in `/cluster/work/climate/dvillanueva/work/CESM/pumas_code/CAM`) and `cam6_4_204` (latest ESCOMP/CAM; CESM tag `cesm3_0_beta10` pins `cam6_4_203`, same microphysics defaults; carries the CAM7 = CESM 3 beta physics; fetched from GitHub, no local checkout). Fetch with `curl -sfL https://raw.githubusercontent.com/ESCOMP/CAM/cam6_4_204/<path>`.
- ☀️🚫 Graupel/rimed ice is **not radiatively active** in either reference: `graupel_in_rad` is a `radiation` namelist switch that defaults to `.false.` in both tags (`namelist_defaults_cam.xml` l.2122 / l.2389), and no CAM7 use case or testmod turns it on. The rel_60 interface therefore matches CESM 3 behaviour as-is.
- 🧊✅ CAM7 (`-phys cam7`) selects `mg3` and sets `micro_mg_do_graupel=.true.`; in this interface graupel stays off unless `micro_mg_do_graupel=.true.` (golden rule). `cnst_add` of `GRAUQM`/`NUMGRA`, `is_convtran1=.true.`, and the `gflx` addition to `LS_FLXPRC`/`LS_FLXSNW` are identical to the upstream host.
- 🎛️📐 CAM7 ships MG3 tuning that is **not** PUMAS default (`micro_mg_vtrmi_factor=1.39`, `micro_mg_accre_enhan_fact=0.731`, `micro_mg_autocon_lwp_exp=2.37`, `micro_mg_implicit_fall`, `micro_mg_accre_sees_auto`, `micro_mg_precip_fall_corr=.true.`). The interface keeps the cam6 defaults; `user_nl_cam.pumas_cesm3like` sets them, together with the CAM7 activation values `microp_aero_wsub_scale=0.9`, `wsubi_scale=1.5`, `wsub_min=0.1`, `wsub_min_asf=0.0`.
- 📡🗑️ Not yet in the interface (only needed if graupel radiation is ever switched on, or for diagnostics): pbuf `DEGRAU`/`ICGRAUWP`/`CLDFGRAU` (filled in `cam6_4_204` `micro_pumas_cam.F90` l.2611-2664), history `LS_REFFGRAU`, `CLDFGRAU`, `UMG`, `FREQG`, `QGSEDTEN`, and the graupel process rates (`PRACGO`, `PSACWGO`, `PGSACWO`, `PGRACSO`, `PRDGO`, `QMULTGO`, `QMULTRGO`, `PSACRO`, `NPRACGO`, `NSCNGO`, `NGRACSO`, `NMULTGO`, `NMULTRGO`, `NPSACWGO`, `MELTGTOT`). Radiative graupel would need overrides of rel_60 `cloud_rad_props.F90` and `radiation.F90` (cam6_3_097 `get_grau_optics_sw`, `grau_cloud_get_rad_props_lw`, `cldfprime = max(cld, cldfgrau)`), behind a `graupel_in_rad` switch defaulting to `.false.`.

## 🧪🔍 Validate

- 🛠️📋 `Buildconf/camconf/Filepath` of the case must start with this directory.
- 🧮✅ A quick standalone compile check of the eleven files, in dependency order, is in `handoff.*.md`.
- 🌍📈 Compare a 3-month run with an MG2 reference in `PRECT`, `CLDLIQ`, `CLDICE`, `TGCLDLWP`, `SWCF`, `LWCF`.

- 🧬⚠️ `create_clone` resets `SRCROOT` to `$W/cesm215/cime/..`. After cloning a PUMAS case, set it back to `$W/cesm215_pumas_worktrees/cesm215_pumas` in `env_case.xml` and `LockedFiles/env_case.xml`, otherwise `build-namelist` rejects the oPUMAS variables.
- ♻️✅ Bit-for-bit check for a new switch: rebuild `pumas_grau`, rerun 1 day, `cmp` the new `cam.r` restart with a saved copy of the old one.

## 📝🧭 Open work

📋🔗 Features that still need interface changes, switches that are wired but untested, and pending runs are tracked in `TODO.md`.

## ⚠️🧭 Caveats

- 🔄🧊 The override `micro_mg_cam.F90` is based on rel_60, not on Diego's `dv_diag` branch. Porting to `dv_diag` means applying `./make_interface.sh diff` to that branch's file.
- 🧮🔁 Two utility modules are initialised with the same constants: `micro_mg_utils` for CAM post-processing and `micro_pumas_utils` inside PUMAS.
- 📦🌐 `module_neural_net` needs the netCDF Fortran module, which the CAM build already provides.
