# 🤝🟣 Handoff plan: oPUMAS.jl, a readable and tested Julia twin of PUMAS (2026-10-02)

🎯📖 Goal: Diego wants to read, understand and modify the PUMAS microphysics in a form that is beautiful and short-lined, with names a general scientist understands.
🧊🏭 The Fortran stays the only code that runs inside CESM; the Julia twin is a reader's edition that also runs offline on saved columns and is tested against the Fortran.
🤖📋 This file is written for the agent that does the rewriting and testing. Work through the milestones in order; each has an acceptance test that must pass before the next starts.

## 0. 🧭📌 Ground rules

- 📘🔗 Read first: `AGENTS.md`, `PHILOSOPHY.md`, `README.md`, `git_work.md`, `tests/README.md`, `translation/README.md`, `code_analysis/PUMAS_stage_map.md`.
- 🌿🟣 Branch: create `dev/julia` from `main_cam_interface` and cherry-pick the commit that holds this file (it was written on `dev/Sppe`). Push only to `origin`.
- 📂🟣 All Julia work lives in `julia/` at the repo root; the Fortran reference driver lives in `tests/column/`. No core `*.F90` file is modified, so `tests/bfb.sh` is not needed for Julia-only commits.
- 🚫🔁 `julia/` is not added to `translation/files.csv`; the name translator does not touch Julia files.
- 🖥️⛔ Nothing runs on `eu-login-*`. Use `herdr-compute` (check `hostname`, `uptime`, `free -g`); it has about 8 CPUs × 2 GB.
- 🟣📦 Julia on Euler: `module load stack/2024-06 julia/1.10.3` (also available: `julia/1.12.4`; confirm the exact prerequisite with `module spider julia/1.10.3`). `~/.julia` already has `NCDatasets` and `CloudMicrophysics`. Use 1.10 (LTS) and commit `Project.toml` and `Manifest.toml`.
- ✅📝 Every commit follows `~/.claude/rules/git-workflow.md` and updates a handoff `handoffs/julia.<shorthash>.md`.
- 🧮🔒 Fidelity rule: the Julia code computes the same arithmetic in the same order as the Fortran (same constants, same limiters, same `min`/`max` clamps, same loop direction in sedimentation). No physics change, no "improvement", no reordering of floating-point operations unless a test proves it is harmless. Suspected Fortran bugs are reproduced and listed in `julia/FINDINGS.md`.

## 1. 🔭✂️ Scope

| In scope (first pass) | Lines | Julia target |
|---|---|---|
| `pumas_kinds.F90` | 16 | `Float64` everywhere, no file |
| `pumas_gamma_function.F90` | 329 | `src/gamma.jl` |
| `micro_pumas_utils.F90` (56 routines) | 3450 | `src/constants.jl`, `src/hydrometeors.jl`, `src/size_distributions.jl`, `src/processes/*.jl` |
| `micro_pumas_v1.F90` (`micro_pumas_init`, `micro_pumas_tend`, sedimentation, `calc_rercld`) | 4997 | `src/parameters.jl`, `src/stages/*.jl`, `src/step.jl` |
| `interface/cesm215/wv_sat_methods.F90` (only the scheme the test case uses) | part | `src/saturation.jl` |
| proc-rates part of `micro_pumas_diags.F90` | part | `src/rates.jl` |

- ⏳🚫 Out of scope until the core is validated: `micro_mg1_0.F90`, TAU bin warm rain (`pumas_stochastic_collect_tau.F90`), the neural-net emulator (`module_neural_net.F90`, `tau_neural_net_quantile.F90`, `ML_fixer_check.F90`), `micro_pumas_ccpp.F90`, `micro_pumas_tempfix.F90`, the MG4/P3 lookup-table routines (`init_lookup_table`, `access_lookup_table*`, `*_mg4`), and the Sppe controls of `dev/Sppe`.
- 🛑📣 For out-of-scope options the Julia code throws a clear error (`error("warm_rain = :tau is not translated yet")`), never a silent fallback.
- 🔀📏 The Fortran has `_line`, `_vect`, `_2D` and `_vec` variants of the same formula. Julia gets one scalar function per formula; arrays use broadcasting or a loop in the stage.

## 2. 🔤📏 Style contract (names and line length)

1. 🌍🔤 Names that cross a function boundary (function names, struct fields, keyword arguments) are full English words from the `CF-CCPP` column of `translation/dictionary.csv` and the AMS Glossary. No unit suffix in Julia names; units go in docstrings and in the dictionary.
2. 🗃️🧩 Structs carry the context, so names do not repeat it: `rates.rain.accretion_by_snow.mass`, not `local_rain_mass_tend_accretion_by_snow`.
3. 🔬🔣 Inside a function, short symbols are allowed (`ρ`, `λ`, `μ`, `N₀`, `qᵣ`, `Δv`) if the unpacking line or the docstring defines each one. Every symbol has a row in `julia/SYMBOLS.md` (symbol, long name, PUMAS name, units).
4. 📏🎨 Maximum line length is 92 columns, enforced by JuliaFormatter (`julia/.JuliaFormatter.toml`, `margin = 92`, `style = "blue"`). An expression that does not fit gets a named intermediate that explains the physics.
5. 🧩📜 One function per process, ideally under 30 lines. Each docstring states: what it computes, units of inputs and outputs, literature reference, and `Fortran twin: <routine> in <file>`.
6. 🎛️🔕 Switches are fields of an immutable `Parameters` struct with the Fortran defaults; namelist spellings (`micro_mg_*`) appear in the docstring of each field so a CESM user can find them.
7. 🚫🧠 No metaprogramming, no clever type hierarchies, no abstract dispatch trees. Plain structs, plain functions, `@kwdef`. Readability beats speed; type stability is kept where it costs nothing.
8. 💬🎨 Comments and docs follow Diego's emoji-first rule; emojis never go inside code or identifiers.
9. 🔤📖 New concepts get a dictionary row first: add a `Julia` column to `translation/dictionary.csv` on `main_cam_interface` (struct path, e.g. `rates.rain.accretion_by_snow.mass`), plus a `Julia` row in `translation/languages.csv` only if `check.py` requires it. `julia/tools/check_names.jl` (milestone M7) verifies that every dictionary `Julia` path exists in the package.

🧪📐 Reference example for the style (shape only; copy the arithmetic from the Fortran twin, not from here):

```julia
"""
Accretion of rain by snow: mass (kg kg⁻¹ s⁻¹) and number (kg⁻¹ s⁻¹) tendencies.
Reference: Ikawa and Saito (1991). Fortran twin: accrete_rain_snow in micro_pumas_utils.F90.
"""
function accretion_of_rain_by_snow(air, rain, snow)
    (; ρ, T) = air
    qᵣ, N₀ᵣ, λᵣ, vᵣ = rain.q, rain.N₀, rain.λ, rain.fall_speed_mass
    qₛ, N₀ₛ, λₛ, vₛ = snow.q, snow.N₀, snow.λ, snow.fall_speed_mass
    ...
    return (; mass, number)
end
```

## 3. 📂🟣 Package layout

```
julia/
  Project.toml  Manifest.toml  .JuliaFormatter.toml
  README.md            how to run, test, read; map Fortran routine -> Julia function
  SYMBOLS.md           symbol table
  FINDINGS.md          oddities found in the Fortran while translating
  TWINS.csv            julia_function, fortran_file, fortran_routine, upstream_commit
  src/
    oPUMAS.jl          module, includes, exports
    constants.jl       physical constants and tunable defaults
    parameters.jl      Parameters struct (switches), built by setup(...) = micro_pumas_init
    hydrometeors.jl    HydrometeorProperties (MGHydrometeorProps), state structs
    gamma.jl  saturation.jl  size_distributions.jl
    rates.jl           ProcessRates (proc_rates_type) and local rate containers
    processes/         one file per family: warm_rain.jl, ice_nucleation.jl, freezing.jl,
                       deposition.jl, aggregation.jl, riming.jl, evaporation.jl, graupel.jl
    stages/            one file per stage of code_analysis/PUMAS_stage_map.md:
                       01_setup.jl ... 15_diagnostics.jl (+ sedimentation.jl)
    step.jl            microphysics_step!(...) = micro_pumas_tend, about 16 calls
  test/
    runtests.jl  test_gamma.jl  test_processes.jl  test_stages.jl  test_step.jl
  tools/
    check_names.jl  check_twins.sh
tests/column/
  dump/                CAM-side dump patch (SourceMods), see M1
  driver.F90           standalone Fortran driver
  build.sh  run.sh  README.md
```

## 4. 🪜🏁 Milestones

### M0 🟣🧰 Bootstrap (small)
- 🌿 Create `dev/julia`, `julia/` skeleton, `Project.toml` (deps: `NCDatasets`, `Test`; dev tool: `JuliaFormatter`), empty module, `runtests.jl` that passes.
- 📘 Add principle 5 to `PHILOSOPHY.md`: "Julia twin: hand-written, verified numerically, never run in CESM, every function names its Fortran twin". Add one line to root `README.md` and one rule to `AGENTS.md` (`julia/` tests before committing Julia changes).
- ✅ Accept: `julia --project=julia -e 'using Pkg; Pkg.test()'` passes on `herdr-compute`.

### M1 🧪📦 Reference data from the Fortran (the critical milestone)
🎯 Produce golden NetCDF files that hold every input and output of `micro_pumas_init` and `micro_pumas_tend`, including all `proc_rates` components.

1. 📤 **Dump inputs in CAM.** In a clone of case `cesm215_pumas.F2000climo.bfb_dev` (use `create_clone --keepexe` as `tests/bfb.sh par` does), add to the SourceMods copy of `interface/cesm215/micro_mg_cam.F90` a debug block, guarded by an environment or namelist switch that is off by default, that writes all arguments of the `micro_pumas_tend` call (after packing to `mgncol` columns) and all arguments of `micro_pumas_init` to a raw binary or NetCDF file for one chunk and a few time steps. Keep this patch in `tests/column/dump/`, not in `interface/cesm215/`, so the golden rule (defaults equal upstream) is untouched.
2. 🔎 **Select columns.** From the dump pick about 200 columns that cover: warm cloud, mixed phase, cirrus, deep convection outflow, polar, clear sky, precipitation without cloud. Check coverage by counting, per `proc_rates` component, how many selected columns have a non-zero value; every in-scope process must be non-zero in at least 5 columns. Add hand-made synthetic columns for processes that stay at zero.
3. 🛠️ **Standalone driver.** `tests/column/driver.F90` reads the dump, calls `micro_pumas_init` and `micro_pumas_tend` from the repo sources, and writes `golden_<label>.nc` with all outputs. `build.sh` compiles `pumas_kinds`, `pumas_gamma_function`, `micro_pumas_sppe` (when on a branch that has it), `wv_sat_methods`, `micro_pumas_utils`, `micro_pumas_diags`, the TAU/ML modules and `micro_pumas_v1` with `gfortran -O0 -fno-fast-math -ffp-contract=off` (module line as in `tests/compile_check.sh`). `-O0` without FMA contraction is what makes near-bitwise agreement with Julia possible.
4. ♻️ **Self-check.** The driver's outputs must equal the outputs dumped by CAM for the same columns (bitwise if the driver is built with the case's flags, otherwise to 1e-12 relative). This proves the driver calls PUMAS the way CAM does.
5. 🎛️ **Switch rows.** Run the driver for every in-scope row of `tests/switches.csv` by changing the init arguments; one golden file per row. Rows for out-of-scope features (`tau`, emulated, MG4, Sppe) are listed as skipped in `tests/column/README.md`.
6. 🔬 **Process-level data.** ✅ Done differently: `tests/column/capture/patch_v1.py` copies `micro_pumas_v1.F90` with a record of every process call inside `micro_pumas_tend` (arguments before and after), so the data are the real calls on the real columns. Base settings go to `golden_processes.nc`, Seifert-Beheng warm rain to `golden_processes_sb2001.nc`. Layout in `tests/column/README.md`.
7. 🪜 **Stage-level data.** Optional but valuable: a copy of `micro_pumas_v1.F90` inside `tests/column/` with dump calls at the 16 stage boundaries of the stage map, so Julia stages can be tested one at a time. The copy is generated by a small patch script so it can be regenerated after upstream merges.
- 🚫📦 Golden files are data: keep them in `tests/column/data/` and add that path to `.gitignore`; `run.sh` regenerates them and prints their absolute paths. Store them on `/cluster/work/climate/dvillanueva/` if large.
- ✅ Accept: self-check (4) passes; coverage table (2) is in `tests/column/README.md`; `run.sh` regenerates everything in under 4 minutes on `herdr-compute`.

### M2 🧱🔢 Foundations
- 🟣 `constants.jl`, `gamma.jl`, `saturation.jl`, `hydrometeors.jl`, `parameters.jl` with `setup(; kwargs...)` mirroring `micro_pumas_init` and `micro_pumas_utils_init` (including derived constants computed there).
- ✅ Gamma: gfortran builds `pumas_gamma` with the intrinsic, which is glibc `tgamma`; `gamma_function` calls the same `tgamma` through `ccall`, so the bits are identical.
- ✅ Accept: gamma and saturation match `golden_processes.nc` to 1e-14 relative; every derived constant of `setup` equals the driver's value bitwise.

### M3 🌧️❄️ Process functions
- 🟣 Translate the in-scope routines of `micro_pumas_utils.F90`, one scalar function each, grouped in `src/processes/`. Order: size distributions, warm rain (`kk2000`, `sb2001v2`, accretion, self-collection), ice nucleation and freezing, deposition/Bergeron, aggregation, riming and secondary ice, evaporation/sublimation, graupel.
- 📋 Add each to `TWINS.csv` when written.
- ✅ Accept per function: max relative error ≤ 1e-13 against `golden_processes.nc`, with exact zeros where the Fortran gives zero (thresholds such as `qsmall` must branch identically).

### M4 🧊🪜 The time step, stage by stage
- 🟣 Design the containers first and get Diego's review on them before filling stages: `ColumnInputs`, `State`, `InCloudValues`, `LocalRates`, `ProcessRates`, `Tendencies`, `Diagnostics`. Arrays are `(column, level)` like the Fortran, level 1 at the top.
- 🟣 Translate stages 1 to 15 of the stage map in order, one file each; `microphysics_step!` calls them in sequence. Stage 8 (conservation limiters) and stage 12 (sedimentation, explicit and implicit, including `implicit_fall`) are the delicate ones: keep loop order and sub-stepping exactly.
- 🧮 The about 280 local work arrays of the Fortran become fields of the containers or locals of a stage; none are module-level globals.
- ✅ Accept per stage (if M1.7 was done): outputs at the stage boundary match to 1e-12 relative.
- ✅ Accept for the step: all outputs and all `proc_rates` of `golden_base.nc` match with max relative error ≤ 1e-10 and absolute error below 1e-20 where the reference is zero; report the worst variable and column in the test log.
- 📌 Status (2026-10-02): done. All 155 fields of `golden_base.nc` are bit-for-bit identical once real powers, exp and log come from libm; per-stage checks (M1.7) were not needed.

### M5 🎛️✅ Switch matrix and conservation
- 🟣 Run `test_step.jl` over every in-scope golden file from M1.5.
- ⚖️ Add physical tests that need no Fortran: total water and energy conservation of one step, no negative mixing ratios, zero tendencies for a clear dry column.
- ✅ Accept: all in-scope switch rows pass the M4 tolerance; conservation holds to rounding.
- 📌 Status (2026-10-02): done. All 16 switch rows with a golden file are bit-for-bit (155 fields, 159 for sb2001; `tau` and `wsub_default` have no golden file). Column water closes to ≤ 3e-18 kg m⁻² s⁻¹ and energy to ≤ 2e-12 W m⁻² (relative 1e-15); negative mixing ratios stay within `qsmall Δt` (a PUMAS limiter threshold, FINDINGS.md); a clear dry column has zero tendencies.

### M6 📖✨ Readability pass (after the numbers are locked)
- 🎨 Run JuliaFormatter; no line over 92 columns (`tools/` check greps for it).
- 🔤 Replace remaining PUMAS-style scratch names (`dum`, `dum1`, `ratio`) by named intermediates; complete `SYMBOLS.md`.
- 📘 Write `julia/README.md`: how to run one column, how to add a process (mirror of the checklist in the stage map), table Fortran routine → Julia function.
- 🧪 Tests must stay green after every renaming commit; renaming and arithmetic changes never share a commit.
- ✅ Accept: Diego reads `step.jl` and two process files and approves the style.
- 📌 Status (2026-10-02): done; Diego approved the style. `s08_conservation.jl` was then split into one function per limiter with a shared `scale!` helper and the rates bound once in `limiter_rates`, still bit-for-bit. Stages bind readable locals with the PUMAS name in a comment (`work.fortran` aliases removed), constants and time step under equation symbols, `dum`-style scratch replaced by named intermediates, JuliaFormatter (sciml, 92 columns) via `tools/format.sh`, `SYMBOLS.md` lists all 349 stage bindings, README has the one-step example, the add-a-process checklist and the twin table. Bit-for-bit after every commit.

### M7 🔗🔄 Staying in sync with upstream
- 🟣 `tools/check_twins.sh <old>..<new>`: lists Fortran routines changed in a commit range (via `git diff -U0` plus routine line ranges) and prints the Julia twins from `TWINS.csv` that need review. Add it to the "Agent routine" in `git_work.md` after an upstream merge.
- 🟣 `tools/check_names.jl`: every `Julia` path in the dictionary exists; every exported function has a `Fortran twin` line.
- ✅ Accept: both tools run clean; a fake change to `accrete_rain_snow` is reported.
- 📌 Status (2026-10-02): done. `check_twins.sh` reports the fake `accrete_rain_snow` change as `rain_collection_by_snow`; `tools/check_names.jl` passes (215 dictionary paths, 35 exported twins). The `Julia` column is metadata in `translate.py`, so no `languages.csv` row; both tools are in the agent routine of `git_work.md`.

### M8 ⏳🧩 Later, only on Diego's request
- 🧊 Sppe controls from `dev/Sppe`; TAU bin warm rain; MG4 lookup tables; an offline single-column time loop with plots; a `CF-CCPP`-names Fortran driver.
- 📌 Status (2026-10-02): 🧊 Sppe done. The 11 controls are `Settings` fields (off-by-default values = upstream), read from the same namelist names as `dev/Sppe`; the golden files come from `dev/Sppe`'s Fortran via `SPPE=1 tests/column/build.sh` without merging that branch. All 18 rows of `tests/column/switches_sppe.csv` are bit-for-bit (5617 tests). 🌧️ TAU done: `src/processes/bin_warm_rain.jl` (four twins of `pumas_stochastic_collect_tau.F90`) and `tau_warm_rain!` in stage 5; the `tau` row is bit-for-bit on all 172 fields, including the 17 TAU process rates (5793 tests). MG4, time loop and CF-CCPP driver remain.

## 5. ⚠️🕳️ Known pitfalls

- 🔢 Fortran `**` with integer exponents, `sign`, `nint`, `int` and integer division differ from naive Julia; use `^` with integer literals, `copysign`, `round(Int, x, RoundNearestTiesAway)` and `div`.
- 🧮 Fortran evaluates `a*b*c` left to right; keep parentheses and order. Do not use `@fastmath`, `muladd` or `fma`.
- 🔁 `where`/masked vector code in `_vect` routines becomes an `if` in the scalar function; check that the `else` branch sets the same zeros.
- 🧊 Module-level saved variables set in `micro_pumas_init` (e.g. `dcs`, `rhmini`, `nccons`, graupel/hail density choice) all become `Parameters` fields; grep `micro_pumas_v1.F90` lines 130–320 for the full list.
- ⬇️ Sedimentation sub-step counts depend on CFL with `int`/`ceiling`; a one-ulp difference can change the count, so compare `nstep` itself in the tests.
- 🌡️ `wv_sat_methods` has several saturation schemes; translate the one selected in the CESM 2.1.5 case (check `wv_sat_init` in the case) and error on the others.
- 🖥️ OpenACC directives (482 lines) are ignored; they carry no physics.
- 📏 `mgncol` packing: the driver and Julia both work on already packed columns; `micro_pumas_get_cols` is not translated.

## 6. 🙋🎛️ Assumptions made for Diego (change here before starting if wrong)

- 🎯 Purpose: readable and runnable offline; Fortran remains production.
- 🔤 Names: `CF-CCPP` vocabulary, struct paths, no unit suffixes, Unicode symbols allowed inside functions, 92 columns.
- 📂 Location: `julia/` on branch `dev/julia` (from `main_cam_interface`), merged back when M5 passes.
- 📏 Agreement target: 1e-10 relative for the full step against a `-O0` gfortran build; reached bit-for-bit in M4/M5 once real powers, exp and log come from libm.

## 7. 📊🧾 What to report at each milestone

- ✅ The acceptance numbers (worst relative error, variable, column), the commands run, absolute paths of files written.
- 🪞 The README reflection required by `AGENT_INSTRUCTIONS.md`.
- 🧠 New entries in `julia/FINDINGS.md`.
