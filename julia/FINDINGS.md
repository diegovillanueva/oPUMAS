# 🧠🔎 Findings from translating PUMAS to Julia

🧾📌 Oddities in the Fortran that the Julia twin keeps on purpose. Each entry names the Fortran place and what the twin does.

## 🧮🔢 Gamma function is the compiler intrinsic, not PUMAS' own series
- 📍 `pumas_gamma_function.F90`: with gfortran `__GFORTRAN__` defines `HAVE_GAMMA_INTRINSICS`, so `pumas_gamma` is the intrinsic `gamma`, which is glibc `tgamma`.
- 🟣 The twin calls `tgamma` from `libm` through `ccall`, not `SpecialFunctions.gamma` and not a port of the series in the file.
- ⚠️ The plan (handoffs/julia.plan.md, M2) assumed the series; it is unused with gfortran builds.

## 📐⚖️ Real powers, exp and log come from libm, as in the Fortran
- 📍 Fortran `x**y` with a real exponent, `exp`, `log` and `log10` call glibc libm; Julia's `^`, `exp` and `log` are Julia's own implementations and can differ by one ulp.
- 🧨 One ulp is enough to break the full step: the implicit sedimentation fluxes cancel to near zero at the lowest levels (for example `lflx` 2.5e-29 against -1.2e-27), so native Julia math gave relative errors up to 50 there.
- 🟣 `src/libm.jl` shadows `^(::Float64, ::Float64)`, `exp`, `log` and `log10` inside the module with `ccall`s to `libm.so.6`. With it every routine and all 155 fields of `golden_base.nc` are bit-for-bit identical.
- 🔢 Integer powers stay Julia's (`x^2` is `x*x`); powers by an integer literal larger than 3 are written out the way gfortran's `powi` multiplies (`x**7` is `x³·x⁴`).

## 🧊📏 Ice lambda bound uses Fortran precedence
- 📍 `micro_pumas_utils.F90`, `micro_pumas_utils_init`: the lower ice lambda bound is written `1._r8/2._r8*dcs`, which Fortran evaluates as `(1/2)*dcs`, not `1/(2*dcs)`.
- 🐛🤔 With dcs = 500e-6 m the bound is 2.5e-4 per metre, so it never binds. MG2 wrote `1._r8/[2._r8*dcs, 1.e-6_r8]`, a bound of 1000 per metre. This looks like an upstream precedence slip.
- 🟣 The twin reproduces `dcs/2` and says so in a comment; a fix belongs upstream, not in the twin.

## 🔢🧊 One constant is folded with higher precision by the compiler
- 📍 `micro_pumas_utils.F90`: `droplet_mass_25um = 4._r8/3._r8*pi*rhow*(25.e-6_r8)**3` is a parameter, so gfortran folds it with exact arithmetic (MPFR) and rounds once.
- 🟣 Julia computes the same expression at run time with rounding after each operation and lands one unit in the last place lower. The twin writes the folded value as a literal; every other folded constant (`aj`, `mi0`, `mmult`, `mi0l_min`, `droplet_mass_40um`) matches bit for bit when computed.

## 🔁💧 CAM passes `allow_sed_supersat` as `remove_supersat`
- 📍 `interface/cesm215/micro_mg_cam.F90` calls `micro_pumas_init(..., allow_sed_supersat, micro_pumas_warm_rain, ...)`, and the dummy argument at that position is `remove_supersat_in`.
- 🐛🤔 A namelist that allows supersaturation therefore asks PUMAS to remove it. Upstream CAM passes the same variable, so the twin keeps the Fortran meaning: `Settings.remove_supersaturation` is read from `allow_sed_supersat` unchanged.

## 🟤🌬️ Graupel sublimation mixes the graupel and snow fall speed exponents
- 📍 `micro_pumas_utils.F90`, `evaporate_sublimate_precip_graupel`: the ventilation term of graupel uses `gamma(5/2 + bg/2)` but divides by `lamg**(5/2 + bs/2)`, the snow exponent. It also uses the snow ventilation coefficients `f1s`, `f2s`.
- 🐛🤔 Dimensionally the λ power should follow the gamma argument, `5/2 + bg/2`. With bg = 0.37 and bs = 0.41 the exponents differ by 0.02, so the rate is off by a factor λ^0.02, about 20 % at λ = 10⁴ m⁻¹.
- 🟣 The twin keeps the snow exponent and marks it in `precipitation_evaporation`.

## 🟤❄️ Rime splintering on graupel has two branches that never run
- 📍 `graupel_rime_splintering`: inside `if (t < tm_3 .and. t > tm_8)` the code still tests `t > tm_3` and `t < tm_8`, both impossible there.
- 🟣 The twin keeps only the two live branches, (T₀ - 3 - T)/2 above -5 °C and (T - T₀ + 8)/3 below, which give the same numbers.

## 🌡️🔢 Some temperature thresholds are hard-coded in kelvin
- 📍 `secondary_ice_production` uses 270.16, 268.16 and 265.16 K, which are -3, -5 and -8 °C from the triple point (273.16 K). `graupel_rime_splintering` uses `tmelt - 3` and so on from 273.15 K. `heterogeneous_rain_freezing` uses 269.15 K.
- 🐛🤔 The two Hallett-Mossop windows of snow and graupel are therefore 0.01 K apart and do not follow a host that changes `tmelt`.
- 🟣 The twin writes the same literals.

## 🧾🌍 The golden inputs hold real NaN values
- 📍 `tests/column/data/golden_inputs.nc` has NaN in 5 points of `tnd_nsnow` and 1 of `re_ice`, and `_FillValue = NaN`, so NCDatasets reads them as `missing`.
- 🟣 `read_inputs` reads the raw values (`.var`) and the step test treats NaN against NaN as equal. The Fortran outputs are reproduced bit for bit, so the NaNs do not reach any output with these settings.

## ➖💧 Conservation limiters skip sink rates below `qsmall`
- 📍 The rain limiter (`micro_pumas_v1.F90` ~l. 2721) only rescales sinks whose summed rate is `>= qsmall`, which compares a rate (kg kg⁻¹ s⁻¹) with a mixing ratio threshold. The snow mass limiter has the same test; cloud liquid, cloud ice and graupel mass are limited whenever the sink exceeds the mixing ratio.
- 🔢 A smaller sink is never limited, so over `Δt` it may remove up to `qsmall Δt` (6e-16 kg kg⁻¹ at 600 s) more than is present. In `golden_inputs.nc` column 156, level 22, rain freezing to ice (`mnuccri` = 2.0e-19 kg kg⁻¹ s⁻¹) meets 1.6e-17 kg kg⁻¹ of rain and leaves −1.07e-16 kg kg⁻¹, in every switch row.
- 🧪 `test_physics.jl` therefore allows negatives down to `qsmall Δt` plus rounding. Column water still closes to ≤ 3e-18 kg m⁻² s⁻¹ and energy to 1e-15 relative: the overdraft is booked in the tendencies, not lost.

## 🧊🪜 Quirks of the time step (stages 11–15 of `micro_pumas_tend`)
- 🧊 The ice fall speed factor `vtrmi_factor` multiplies `vtrmi` and then again the blended small-ice speed when the ice fraction is below 1, so blended points get the factor twice.
- 🪂 The precipitation fall correction starts at `k > 2`, so the second level never takes the speed of the level above.
- 🟤 In the effective radius stage `dumg` and `dumng` are grid means, not divided by `precip_frac` like rain and snow; the graupel number check then multiplies by `precip_frac` again.
- 🟤 The 10 cm reflectivity of graupel uses `n0g` from the process stage, which stage 15 never updates.
- 🌧️ The rain freezing after sedimentation sizes the rain from grid-mean `dumr`, not in-precipitation values.
- 🟤 `do_hail .or. do_graupel .and. lamg > 0` binds as `do_hail .or. (do_graupel .and. lamg > 0)`; with hail on, graupel reflectivity is added even when `lamg` is 0.
- 🔢 The explicit sedimentation step `rnstep = 1._r8/real(nstep)` divides by a default-kind (single precision) real; exact for every realistic `nstep`, kept as `Float64(Float32(nstep))`.
- 🚫 `nrcons` and `nscons` (constant rain and snow number) are not settings of oPUMAS, so their branches are left out; they are false in all CAM namelists.
