# 🟣📖 oPUMAS.jl: the readable Julia twin of PUMAS

🌧️🧊 A hand-written Julia version of the PUMAS cloud microphysics (`micro_pumas_v1.F90`), with names a general scientist can read.
⚖️🧪 It is verified against the Fortran with a standalone column driver ([`../tests/column/`](../tests/column/)); CESM keeps running the Fortran.

## 🚀 Run

```bash
julia/tools/julia.sh -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

- 🧩🔧 Always use `julia/tools/julia.sh`: it calls Julia 1.10.3 by absolute path (`module load julia` does not reach the PATH here) and clears Euler's `LD_LIBRARY_PATH`, whose gcc-8 `libgfortran` breaks `NetCDF_jll`.
- 🖥️⛔ Run on `herdr-compute`, not on login nodes.
- 📦🧷 The tests read the golden files of `tests/column/data/` (git-ignored). Make them once with `tests/column/run.sh` (under a minute), or point `OPUMAS_GOLDEN` at another directory holding `golden_inputs.nc`, `golden_processes.nc`, `golden_processes_sb2001.nc` and `settings_<label>.nml` plus `golden_<label>.nc` for the switch rows of `../tests/switches.csv` (rows without a golden file, `tau` and `wsub_default`, are skipped).
- 🔬📋 The tests print, per switch row, how many step fields are bit-for-bit identical (all 155, 159 for `sb2001`, since M5) and the column water and energy budget residuals, and end with the worst relative error per Fortran routine (process twins held to 1e-13).
- ⏱️ `Pkg.test()` takes about 90 s on herdr-compute, mostly compilation; one warm step of the 209 columns takes about 30 ms.
- 📐🧮 Bit-for-bit needs libm: `src/libm.jl` routes real powers, `exp`, `log` and `log10` to glibc as gfortran does (FINDINGS.md).

## 📂 Layout

- 📦 `src/`: the package (`oPUMAS.jl` includes everything; `libm.jl` first). Foundations (constants, saturation, hydrometeors, parameters, namelist, size distributions) at the top, one file per process family in `src/processes/`: `warm_rain`, `aggregation` (ice and snow), `riming` (snow collecting droplets and rain), `freezing`, `deposition` (vapor growth, Bergeron, evaporation), `graupel`, `effective_radius`.
- 🪜 The time step: `names.jl` (PUMAS name → container path tables, `output_paths()`), `workspace.jl` (`ColumnInputs`, `read_inputs`, `allocate`, and the PUMAS-name aliases `work.fortran`), `stages/s01_…s15_*.jl` (one file per stage of `micro_pumas_tend`, short names in `stages/common.jl`) and `step.jl` (`microphysics_step!`).
- ▶️ One step: `work = allocate(read_inputs(path)); microphysics_step!(work, setup(read_settings(nml)))`.
- 🧪 `test/`: tests against the golden files from `tests/column/`. `processes.jl` replays captured Fortran calls; one `test_<topic>.jl` per topic (`test_ice_processes.jl` covers every ice, snow, graupel and vapor process; `test_step.jl` runs the full step against every `golden_<label>.nc`; `test_physics.jl` checks water and energy budgets, negative mixing ratios and a clear dry column without the Fortran).
- 📋 `TWINS.csv`: Julia function, Fortran file and routine, upstream commit translated.
- 🔣 `SYMBOLS.md`: every short symbol used inside functions.
- 🐛 `FINDINGS.md`: oddities of the Fortran met while translating.
- 🛠️ `tools/`: the `julia.sh` wrapper and `check_twins.sh <old>..<new>`, which lists the twins whose Fortran routines changed between two commits.
- 📐🧱 `DESIGN.md`: the M4 container design (draft answers adopted; the stages still use PUMAS names until M6).
- 🗺️ Plan and milestones: [`../handoffs/julia.plan.md`](../handoffs/julia.plan.md).
