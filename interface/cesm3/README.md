# 🌍🧊 oPUMAS Sppe controls for CESM 3 / CAM7

🎯📦 This directory carries the CAM-side part of the Sppe controls (namelist group `micro_dv_nl`, from
[Sppe_revision_public](https://github.com/diegovillanueva/Sppe_revision_public)) for CESM 3 / CAM7.
🧊🔌 The PUMAS-side controls live in oPUMAS itself (`micro_pumas_sppe.F90`, hooks in `micro_pumas_v1.F90` and
`micro_pumas_utils.F90`), so CESM 3 only needs these patches plus oPUMAS as its PUMAS source.

⚠️🧪 **Status: checked, not compiled or run.** There is no CESM 3 checkout on Euler. The patches are checked
with `git apply --check` against CAM tag `cam6_4_204`; the tested implementation is `../cesm215`.

## 📂 Contents

| File | Purpose |
|---|---|
| `patches/micro_dv_nml.F90.patch` | 🆕 new module: reads `micro_dv_nl`, broadcasts it (same as `../cesm215/micro_dv_nml.F90` minus `use mpishorthand`) |
| `patches/microp_aero.F90.patch` | 📥 `call micro_dv_readnl(nlfile)` in `microp_aero_readnl` |
| `patches/micro_pumas_cam.F90.patch` | 📞 `micro_pumas_sppe_set(...)` right after `micro_pumas_init` in `micro_pumas_cam_init` |
| `patches/hetfrz_classnuc.F90.patch` | ❄️ `depr_point_frz` and `limfacdu_niemand_continuous`, rewritten for `hetfrz_classnuc_calc_rates` |
| `patches/hetfrz_classnuc_cam.F90.patch` | 🎛️ `hclas_{bc,du}_{imm,dep,cnt}` factors on the summed rates |
| `patches/clubb_intr.F90.patch`, `patches/macrop_driver.F90.patch` | ☁️ `detrainment_ramp_liq` in the detrainment phase ramp |
| `patches/namelist_definition.xml.patch` | 📜 the 20 `micro_dv_nl` entries |
| `user_nl_cam.sppe` | 🧊 the full Sppe setting for CAM7 |
| `make_interface.sh` | 🧪 `check` (dry run against the tag) and `apply <CAMROOT>` (idempotent) |

## ▶️ How to use

```bash
interface/cesm3/make_interface.sh check                       # 🧪 needs network (raw.githubusercontent.com)
interface/cesm3/make_interface.sh apply <cesm3>/components/cam # 🩹 in a CAM7 checkout at cam6_4_204
```

🔗🧊 Then point CAM's PUMAS submodule (`src/physics/pumas-frozen` in `.gitmodules`, `fxtag =
pumas_cam-release_v1.17_rename`) at oPUMAS `dev/Sppe`, so `micro_pumas_sppe.F90` and the core hooks are compiled.

## 🗺️ Differences from CESM 2.1.5

| Control | CAM7 (`cam6_4_204`) |
|---|---|
| `limfacdu` | 🔁 read but unused: CAM7 already has `hetfrz_dust_scalfac` (group `hetfrz_classnuc_nl`, default 0.05). It also caps dust deposition and contact freezing, not only immersion as `limfacdu` does in rel_60. |
| `depr_point_frz = .false.` | ❄️ skips the `aw*supersatice > 1` gate in `hetfrz_classnuc_calc_rates` and keeps the unmodified germ radius, as in rel_60 |
| `limfacdu_niemand_continuous` | 🧮 applies to the aerosol types with `pdf_imm` (dust) and uses each type's `mradius`, not the fixed rel_60 dust-mode radii |
| `detrainment_ramp_liq` | ☁️ in `clubb_intr` the ramp lower bound is `clubb_detphase_lowtemp`, not a fixed 238.15 K; the switch still sets `dum1 = 0` inside the ramp |

## 🚧 Open

- ⚠️🔢 `cam6_4_204` pins PUMAS `v1.17_rename`; oPUMAS is v1.41+. Whether `micro_pumas_cam.F90` of that tag compiles
  against oPUMAS (argument lists of `micro_pumas_init`/`micro_pumas_tend`) is not checked.
- 🧪 No bit-for-bit or switch test exists for CESM 3.
