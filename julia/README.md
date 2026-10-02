# 🟣📖 oPUMAS.jl: the readable Julia twin of PUMAS

🌧️🧊 A hand-written Julia version of the PUMAS cloud microphysics (`micro_pumas_v1.F90`), with names a general scientist can read.
⚖️🧪 It is verified against the Fortran with a standalone column driver ([`../tests/column/`](../tests/column/)); CESM keeps running the Fortran.

## 🚀 Run

```bash
julia/tools/julia.sh -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

- 🧩🔧 Always use `julia/tools/julia.sh`: it calls Julia 1.10.3 by absolute path (`module load julia` does not reach the PATH here) and clears Euler's `LD_LIBRARY_PATH`, whose gcc-8 `libgfortran` breaks `NetCDF_jll`.
- 🖥️⛔ Run on `herdr-compute`, not on login nodes.

## 📂 Layout

- 📦 `src/`: the package (`oPUMAS.jl` includes everything).
- 🧪 `test/`: tests against the golden files from `tests/column/`.
- 🛠️ `tools/`: the wrapper and the name/twin checks.
- 🗺️ Plan and milestones: [`../handoffs/julia.plan.md`](../handoffs/julia.plan.md).
