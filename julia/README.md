# 🟣📖 oPUMAS.jl: the readable Julia twin of PUMAS

🌧️🧊 A hand-written Julia version of the PUMAS cloud microphysics (`micro_pumas_v1.F90`), with names a general scientist can read.
⚖️🧪 It is verified against the Fortran with a standalone column driver ([`../tests/column/`](../tests/column/)); CESM keeps running the Fortran.

## 🚀 Run

```bash
julia/tools/julia.sh -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

- 🧩🔧 Always use `julia/tools/julia.sh`: it calls Julia 1.10.3 by absolute path (`module load julia` does not reach the PATH here) and clears Euler's `LD_LIBRARY_PATH`, whose gcc-8 `libgfortran` breaks `NetCDF_jll`.
- 🖥️⛔ Run on `herdr-compute`, not on login nodes.
- 📦🧷 The tests read the golden files of `tests/column/data/` (git-ignored). Make them once with `tests/column/run.sh` (under a minute), or point `OPUMAS_GOLDEN` at another directory holding `settings_base.nml`, `golden_processes.nc` and `golden_processes_sb2001.nc`.
- 🔬📋 The tests end with the worst relative error per Fortran routine; process twins are held to 1e-13 (1e-12 where Fortran's real powers differ from Julia's).

## 📂 Layout

- 📦 `src/`: the package (`oPUMAS.jl` includes everything). Foundations (constants, saturation, hydrometeors, parameters, namelist, size distributions) at the top, one file per process family in `src/processes/`: `warm_rain`, `aggregation` (ice and snow), `riming` (snow collecting droplets and rain), `freezing`, `deposition` (vapor growth, Bergeron, evaporation), `graupel`, `effective_radius`.
- 🧪 `test/`: tests against the golden files from `tests/column/`. `processes.jl` replays captured Fortran calls; one `test_<topic>.jl` per topic (`test_ice_processes.jl` covers every ice, snow, graupel and vapor process).
- 📋 `TWINS.csv`: Julia function, Fortran file and routine, upstream commit translated.
- 🔣 `SYMBOLS.md`: every short symbol used inside functions.
- 🐛 `FINDINGS.md`: oddities of the Fortran met while translating.
- 🛠️ `tools/`: the `julia.sh` wrapper and `check_twins.sh <old>..<new>`, which lists the twins whose Fortran routines changed between two commits.
- 📐🧱 `DESIGN.md`: the M4 container design (draft, waiting for Diego's review).
- 🗺️ Plan and milestones: [`../handoffs/julia.plan.md`](../handoffs/julia.plan.md).
