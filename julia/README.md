# 🟣📖 oPUMAS.jl: the readable Julia twin of PUMAS

🌧️🧊 A hand-written Julia version of the PUMAS cloud microphysics (`micro_pumas_v1.F90`), with names a general scientist can read.
⚖️🧪 It is verified against the Fortran with a standalone column driver ([`../tests/column/`](../tests/column/)); CESM keeps running the Fortran.

## 🚀 Run

```bash
julia/tools/julia.sh -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

- 🧩🔧 Always use `julia/tools/julia.sh`: it calls Julia 1.10.3 by absolute path (`module load julia` does not reach the PATH here) and clears Euler's `LD_LIBRARY_PATH`, whose gcc-8 `libgfortran` breaks `NetCDF_jll`.
- 🖥️⛔ Run on `herdr-compute`, not on login nodes.
- 📦🧷 The tests read the golden files of `tests/column/data/` (git-ignored). Make them once with `tests/column/run.sh` (under a minute), or point `OPUMAS_GOLDEN` at another directory holding `golden_inputs.nc`, `golden_processes.nc`, `golden_processes_sb2001.nc` and `settings_<label>.nml` plus `golden_<label>.nc` for the switch rows of `../tests/switches.csv` and `../tests/column/switches_sppe.csv` (rows without a golden file, such as `wsub_default`, are skipped; `settings_tau.nml` names the kernel file `KBARF_tau_kernel.dat` by absolute path).
- 🔬📋 The tests print, per switch row, how many step fields are bit-for-bit identical (all 155, 159 for `sb2001`, since M5) and the column water and energy budget residuals, and end with the worst relative error per Fortran routine (process twins held to 1e-13).
- ⏱️ `Pkg.test()` takes about 60 s on herdr-compute, mostly compilation; one warm step of the 209 columns takes about 30 ms.
- 🎨📏 `julia/tools/format.sh` formats `src/` and `test/` (JuliaFormatter, `.JuliaFormatter.toml`, sciml style, 92 columns); `--check` only reports. Run the check before committing.
- 📐🧮 Bit-for-bit needs libm: `src/libm.jl` routes real powers, `exp`, `log` and `log10` to glibc as gfortran does (FINDINGS.md).

## 📂 Layout

- 📦 `src/`: the package (`oPUMAS.jl` includes everything; `libm.jl` first). Foundations (constants, saturation, hydrometeors, parameters, namelist, size distributions) at the top, one file per process family in `src/processes/`: `warm_rain`, `aggregation` (ice and snow), `riming` (snow collecting droplets and rain), `freezing`, `deposition` (vapor growth, Bergeron, evaporation), `graupel`, `effective_radius`.
- 🪜 The time step: `names.jl` (PUMAS name → container path tables, `output_paths()`), `workspace.jl` (`ColumnInputs`, `read_inputs`, `allocate`), `stages/s01_…s15_*.jl` (one file per stage of `micro_pumas_tend`; each binds its fields to readable locals with the PUMAS name in a comment, `qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc`; constants and densities in `stages/common.jl`) and `step.jl` (`microphysics_step!`).
- 🧪 `test/`: tests against the golden files from `tests/column/`. `processes.jl` replays captured Fortran calls; one `test_<topic>.jl` per topic (`test_ice_processes.jl` covers every ice, snow, graupel and vapor process; `test_step.jl` runs the full step against every `golden_<label>.nc`; `test_physics.jl` checks water and energy budgets, negative mixing ratios and a clear dry column without the Fortran).
- 📋 `TWINS.csv`: Julia function, Fortran file and routine, upstream commit translated.
- 🔣 `SYMBOLS.md`: every short symbol used inside functions, every stage local with its container path and PUMAS name.
- 🐛 `FINDINGS.md`: oddities of the Fortran met while translating.
- 🛠️ `tools/`: the `julia.sh` wrapper, `format.sh`, `check_twins.sh <old>..<new>`, which lists the twins whose Fortran routines changed between two commits, and `check_names.jl` (`julia/tools/julia.sh --project=julia julia/tools/check_names.jl`), which checks that every `Julia` path of `translation/dictionary.csv` exists and every exported twin names its Fortran routine in a `Fortran twin:` docstring line.
- 📐🧱 `DESIGN.md`: the container design of the time step and how it was built.
- 🗺️ Plan and milestones: [`../handoffs/julia.plan.md`](../handoffs/julia.plan.md).

## ▶️🌧️ Run one step on your own columns

```julia
# from the repo root: julia/tools/julia.sh --project=julia
using oPUMAS
work = allocate(read_inputs("tests/column/data/golden_inputs.nc"))  # (column, level) arrays
parameters = setup(read_settings("tests/column/data/settings_base.nml"))
microphysics_step!(work, parameters)
work.tendencies.heating           # tlat, W kg⁻¹
work.diagnostics.surface          # prect, preci as total, ice (one value per column)
field_at(work, Dict(output_paths())[:qrout])  # any output by its PUMAS name
```

- 📥 `read_inputs` reads a file in the layout of `golden_inputs.nc` (`tests/column/run.sh` writes one from the CAM column dump); for other data call `ColumnInputs(; time_step, t, q, qcn, ...)` with (column, level) matrices, level 1 at the top.
- 🎛️ `read_settings` reads a `&pumas_settings` namelist as `tests/column/run.sh` writes it (`name = value` per line, CAM namelist names); entries left out keep the `Settings` defaults, unknown entries are an error.
- 🔁 The workspace is allocated once and overwritten by every call, so a time loop calls `microphysics_step!` again after updating `work.inputs` in place.

## ➕🧩 Add a process

🪞 Mirror of the Fortran checklist in [`../code_analysis/PUMAS_stage_map.md`](../code_analysis/PUMAS_stage_map.md); do the Fortran first, then the twin.

1. 🔤 Add the rate's names to `translation/dictionary.csv` on `main_cam_interface` (PUMAS name, CF-CCPP name, `Julia` path such as `rates.rain.accretion_by_snow.mass`).
2. 🧮 Write the rate as a scalar function in `src/processes/<family>.jl` with a `Fortran twin:` line in its docstring, export it, and add a row to `TWINS.csv`.
3. 🗂️ Add its path to `LOCAL_RATE_PATHS` in `src/names.jl` (and to `PROCESS_RATE_PATHS` if it goes to history); `allocate` builds the field.
4. 📞 Call it in stage 5, 6 or 7 behind its `Settings` switch, off by default, binding a local with the PUMAS name in a comment.
5. ⚖️ Add it to the limiter of every species it removes (`s08_conservation.jl`).
6. ➕ Add it to the tendencies and latent heating (`s09_tendencies.jl`).
7. 🎚️ Read the switch in `src/namelist.jl`, add a row to `../tests/switches.csv`, regenerate the golden files and run `Pkg.test()`: the new row must be bit-for-bit, `test_physics.jl` must still close the budgets.
8. 🔣 Add new locals to `SYMBOLS.md` and run `tools/format.sh --check`.

## 🔗📋 Fortran routine → Julia function

| Fortran file | Fortran routine | Julia function |
|---|---|---|
| `libm` | `tgamma (gfortran gamma intrinsic)` | `gamma_function` |
| `wv_sat_methods.F90` | `wv_sat_qsat_water` | `saturation_over_water` |
| `wv_sat_methods.F90` | `wv_sat_qsat_ice` | `saturation_over_ice` |
| `micro_pumas_utils.F90` | `NewMGHydrometeorProps` | `hydrometeor_properties` |
| `micro_pumas_v1.F90` | `micro_pumas_init` | `setup` |
| `micro_pumas_utils.F90` | `micro_pumas_utils_init` | `setup` |
| `micro_pumas_utils.F90` | `rising_factorial` | `rising_factorial` |
| `micro_pumas_utils.F90` | `var_coef` | `variance_coefficient` |
| `micro_pumas_utils.F90` | `calc_ab` | `latent_heat_factor` |
| `micro_pumas_utils.F90` | `size_dist_param_basic` | `size_distribution` |
| `micro_pumas_utils.F90` | `size_dist_param_liq` | `droplet_size_distribution` |
| `micro_pumas_utils.F90` | `avg_diameter_vec` | `mean_diameter` |
| `micro_pumas_utils.F90` | `kk2000_liq_autoconversion` | `kk2000_autoconversion` |
| `micro_pumas_utils.F90` | `sb2001v2_liq_autoconversion` | `sb2001_autoconversion` |
| `micro_pumas_utils.F90` | `sb2001v2_accre_cld_water_rain` | `sb2001_accretion` |
| `micro_pumas_utils.F90` | `accrete_cloud_water_rain` | `accretion` |
| `micro_pumas_utils.F90` | `self_collection_rain` | `rain_self_collection` |
| `micro_pumas_utils.F90` | `ice_autoconversion` | `ice_to_snow_autoconversion` |
| `micro_pumas_utils.F90` | `snow_self_aggregation` | `snow_self_aggregation` |
| `micro_pumas_utils.F90` | `accrete_cloud_ice_snow` | `ice_collection_by_snow` |
| `micro_pumas_utils.F90` | `accrete_cloud_water_snow` | `droplet_collection_by_snow` |
| `micro_pumas_utils.F90` | `secondary_ice_production` | `rime_splintering_by_snow` |
| `micro_pumas_utils.F90` | `accrete_rain_snow` | `rain_collection_by_snow` |
| `micro_pumas_utils.F90` | `heterogeneous_rain_freezing` | `rain_immersion_freezing` |
| `micro_pumas_utils.F90` | `bergeron_process_snow` | `bergeron_to_snow` |
| `micro_pumas_utils.F90` | `vapor_deposition_onto_snow` | `snow_deposition` |
| `micro_pumas_utils.F90` | `ice_deposition_sublimation` | `ice_deposition` |
| `micro_pumas_utils.F90` | `evaporate_sublimate_precip_graupel` | `precipitation_evaporation` |
| `micro_pumas_utils.F90` | `graupel_collecting_snow` | `snow_collection_by_rain` |
| `micro_pumas_utils.F90` | `graupel_collecting_cld_water` | `droplet_collection_by_graupel` |
| `micro_pumas_utils.F90` | `graupel_riming_liquid_snow` | `rimed_snow_to_graupel` |
| `micro_pumas_utils.F90` | `graupel_collecting_rain` | `rain_collection_by_graupel` |
| `micro_pumas_utils.F90` | `graupel_rain_riming_snow` | `rain_snow_collisions_to_graupel` |
| `micro_pumas_utils.F90` | `graupel_rime_splintering` | `rime_splintering_by_graupel` |
| `micro_pumas_v1.F90` | `calc_rercld` | `add_rain_cloud_effective_radius` |
| `micro_pumas_v1.F90` | `micro_pumas_tend` | `microphysics_step!` |
| `micro_pumas_v1.F90` | `implicit_fall` | `implicit_fall!` |
| `micro_pumas_v1.F90` | `Sedimentation_implicit` | `sediment_implicit!` |
| `micro_pumas_v1.F90` | `Sedimentation` | `sediment_explicit!` |
