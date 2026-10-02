# 🤝🟣 Handoff: oPUMAS.jl after M8 (2026-10-02)

- 🔖 **Commit:** `ecaf4cb` ✨ feat: M8 column driver built from the CF-CCPP translation (branch `dev/julia`, merged into `main_cam_interface`, synced to `CF-CCPP`, all pushed to `origin`).
- ✅ **What was done:**
  - ⚖️ `6e22545`: `s08_conservation.jl` rewritten as one function per limiter, after Diego's review; issue #4 on diegovillanueva/oPUMAS.
  - 🧊 `eff1d16`: the 11 Sppe controls of `dev/Sppe` as `Settings` fields; 18 rows of `tests/column/switches_sppe.csv` bit-for-bit, built from `dev/Sppe` with `SPPE=1 build.sh` and not merged.
  - 🌧️ `4a25023`: TAU bin warm rain (`src/processes/bin_warm_rain.jl`, `tau_warm_rain!` in stage 5); the `tau` row is bit-for-bit on 172 fields.
  - ⏱️ `29c569a`: `time_loop!`/`advance_state!` and the twin `nsteps` argument of `tests/column/driver.F90`; 48 steps bit-for-bit (`golden_loop.nc`); `julia/tools/column_loop.jl` and `plot_column_loop.py` record and plot one day.
  - 🔤 `ecaf4cb`: `tests/column/cf_ccpp_check.sh` builds the driver from the CF-CCPP translation, and all 17 switch rows plus the loop are bit-for-bit. `compare.py` now treats equal infinities as equal, and `check.py` scans `tests/column/*.F90`.
  - 📊🚫 MG4 was dropped because its routines have no call site (`julia/FINDINGS.md`).
- 📌 **Current state:**
  - 🏁 M0–M8 are done. `Pkg.test()` passes 5949 of 5949 tests in about 1.5 min on herdr-compute.
  - 🔧 The capture build (`CAPTURE=1`) is still PUMAS-names only; it does not build on `CF-CCPP`.
- ⏭️ **Next steps:**
  1. 🙋 Diego decides whether the time loop should take time-varying host forcing, e.g. a sequence of CAM dumps, instead of fixed forcing.
  2. 🌍 Optional: dump more CAM steps (e.g. 24 consecutive steps of one chunk) to compare the offline loop with CAM's own trajectory.
- 🧭 **Key decisions:**
  - 🔁 The time loop uses forward Euler, and the Fortran driver applies the same update, so both sides stay bit-for-bit. Negative values are zeroed with a `where`/`<` test, so NaN is left as NaN.
  - 🔤 The CF-CCPP check translates the working tree, not the branch, so it runs before a commit. `sync.sh` produces the same text on `CF-CCPP`.
