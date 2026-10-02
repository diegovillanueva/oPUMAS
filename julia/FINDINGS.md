# 🧠🔎 Findings from translating PUMAS to Julia

🧾📌 Oddities in the Fortran that the Julia twin keeps on purpose. Each entry names the Fortran place and what the twin does.

## 🧮🔢 Gamma function is the compiler intrinsic, not PUMAS' own series
- 📍 `pumas_gamma_function.F90`: with gfortran `__GFORTRAN__` defines `HAVE_GAMMA_INTRINSICS`, so `pumas_gamma` is the intrinsic `gamma`, which is glibc `tgamma`.
- 🟣 The twin calls `tgamma` from `libm` through `ccall`, not `SpecialFunctions.gamma` and not a port of the series in the file.
- ⚠️ The plan (handoffs/julia.plan.md, M2) assumed the series; it is unused with gfortran builds.

## 📐⚖️ Real powers, exp and log are libm in Fortran and native in Julia
- 📍 Fortran `x**y` with a real exponent calls libm `pow`; Julia `^`, `exp` and `log` are Julia's own implementations.
- 🟣 Results may differ by one ulp; the tests therefore use a relative tolerance (1e-12 for single routines, 1e-10 for a full step) instead of bit-for-bit.
- 🔢 `x**2` with an integer literal is `x*x` in both languages; the twin keeps integer literals where the Fortran has them.

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
