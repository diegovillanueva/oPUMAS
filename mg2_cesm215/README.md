# 🧊🛠️ MG2 of CESM 2.1.5 with the dust-ice limiter patches

🎯🧬 This directory holds the two microphysics files of CESM 2.1.5 (CAM tag `cam_cesm2_1_rel_60`, `a03b84b7`, MG2 before it was renamed PUMAS) with the minimal patch set of [mg2_limiter_update](https://github.com/diegovillanueva/mg2_limiter_update) (branch `minimal-patches`).
🧾🌿 The first commit of branch `mg2_limiter_update` adds the stock rel_60 files, the second applies the patches, so `git diff` between them shows the change.

| File | Content |
|---|---|
| `micro_mg2_0.F90` | 🔢 ice-number cap (`enable_nimax`, `naai_het_also_in_mpc`), Bergeron, autoconversion and riming factors, `detrainment_ramp_liq` |
| `micro_mg_utils.F90` | ❄️ contact freezing, Hallett–Mossop (`onoff_sip`, `factor_HM`), rain freezing switches |
| `patches/*.patch` | 🩹 the same two diffs vs `cam_cesm2_1_rel_60` as in mg2_limiter_update |

✅🔁 The patched files are byte-identical to `src/physics/cam/` of a rel_60 tree with all 11 mg2_limiter_update patches applied (`cmp`).

## 🧩📋 CAM-side patches that stay in mg2_limiter_update

🔗🧊 Both files `use micro_dv_nml`, so they compile only together with the CAM-side patches:

| Patch | Role |
|---|---|
| `micro_dv_nml.F90` (new) | 🎛️ namelist group `micro_dv_nl` (23 variables), uses CAM utilities |
| `namelist_definition.xml` | 📝 definitions of `micro_dv_nl` |
| `microp_aero.F90` | 📞 calls the `micro_dv_nl` reader |
| `hetfrz_classnuc.F90`, `hetfrz_classnuc_cam.F90` | 🧊 dust INP limiter, Niemand schemes, freezing-rate factors |
| `nucleate_ice.F90` | 🛡️ `log(0)` guard |
| `micro_mg_cam.F90` | 📊 history diagnostics |
| `clubb_intr.F90`, `macrop_driver.F90` | 🌡️ `detrainment_ramp_liq` |

## 🚀⚙️ Use

```bash
# 🩹 apply the CAM-side patches of mg2_limiter_update, then drop these two files into the case
cp micro_mg2_0.F90 micro_mg_utils.F90 /path/to/case/SourceMods/src.cam/
```
