# 🧊🌍 oPUMAS in MESSy (EMAC and ICON)

🎯 How the MESSy ports call PUMAS, what the host must provide, and the pitfalls found while porting (2026-10-02/03). The ports live in the private repo `github.com/diegovillanueva/messy_pumas`; the PUMAS core files are linked unchanged from this repository (submodule `external/oPUMAS`).

## 🔗 Files and how they are used

- 🧱 Core, linked into `messy/smcl/` under `.f90` names (MESSy's CMake globs `smcl/*.f90` only): `pumas_kinds`, `pumas_gamma_function`, `micro_pumas_utils`, `micro_pumas_diags`, `pumas_stochastic_collect_tau`, `module_neural_net`, `tau_neural_net_quantile`, `ML_fixer_check`, `micro_pumas_v1`, plus `interface/cesm215/wv_sat_methods.F90`.
- 🧩 Glue (MESSy-specific, in `messy_pumas`): EMAC `src/messy_cloud_pumas.f90` (scheme `cloud_param = 7` inside the CLOUD submodel); ICON `messy_pumas.f90` + `messy_pumas_si.f90` (own submodel PUMAS).
- 🚫 Not linked: `micro_pumas_ccpp`, `micro_pumas_tempfix`, `micro_mg1_0` (CCPP or CAM clashes); `dgemm` is only needed for `warm_rain = 'emulated'`.

## 📥 What the host must provide to `micro_pumas_tend`

| Input | EMAC glue | Note |
|---|---|---|
| `t, q, qc, qi` | state at end of step (m1 + tendency*dt) | leapfrog: `deltat = ztmst` (2 dt) |
| `nc` (1/kg) | ECHAM `acdnc`/rho * liquid fraction, or constant | droplet number is prescribed; `nctend` ignored |
| `ni, qr, qs, nr, ns` | advected MESSy tracers `NI/QR/QS/NR/NS_pumas` | rain/snow sediment inside PUMAS, `I_sedi` off |
| `cldn, liqcldf, icecldf` | liquid fraction from a uniform-PDF condensation step (width from ECHAM `rhc`), ice fraction = host `aclc` with a floor where ice exists | PUMAS does no condensation and no cover |
| `naai` (1/kg/s) | Meyers (1992) at the actual ice supersaturation capped at water saturation, only where supersaturated and `icecldf > 0`, number capped at `nin_max` | the core's `micro_mg_icenuc_use_meyers` flag does nothing |
| `npccn, rndst, nacon, frz*` | 0 | no activation tendency, no contact freezing (aerosol constant) |
| `pint, pdel, p` | ECHAM interfaces, **top interface forced positive** | see pitfall 1 |
| `relvar, accre_enhan, qsatfac` | constants 2, 1, 1 | CAM values without CLUBB |

## ⚠️ Pitfalls

1. 🧊 **0 Pa model top kills sedimentation.** With `micro_mg_implicit_fall` the layer heights come from `log(pint(k+1)/pint(k))`; a 0 Pa top interface (ECHAM5) makes the top layer infinitely thick, every flux is exactly 0 and rain/snow mass is lost silently. Guard: `pint(:,1) = MAX(pint(:,1), 0.5*p(:,1))`. ICON's top is about 40 hPa and is safe.
2. ❄️ **Meyers at water saturation nucleates far too much ice** (about 2e5 m-3 at -37 C) in clear supersaturated air; use the actual S_i, require supersaturation, cap at 1e5 m-3 (EMAC IWP 25.5 to 19.5 g/m2, reference 20.5).
3. 📏 **`effi` blows up** (1/lami) where `ni -> 0`; clip effective radii before output (ICON glue: `reffc` 1-50 um, `reffi` 5-300 um).
4. 🔤 `warm_rain` and `micro_mg_precip_frac_method` must be `character(len=16)` actuals (ifort error #7938 otherwise).
5. 🧮 `dcs` (ice-to-snow autoconversion size) is host dependent: ICON needs 350e-6 to match its two-moment reference, EMAC is right at 500e-6 after the nucleation fix.

## 🧪 Validation

- 🌍 EMAC: `messy_pumas/scripts/make_pumas_min_runscripts.sh` and `compare_pumas_min.py` (3-day T42L31 runs, water budget, supercooled liquid fraction); results in `messy_pumas/figures/pumas_min_3day.md`.
- 🌐 ICON: `messy_pumas/figures/pumas_icon_3day.md`.
- ⚡ Compile check without a MESSy build: `messy_pumas/scripts/compile_check.sh`.
