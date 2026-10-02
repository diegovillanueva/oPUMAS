# 🤝🟣 Handoff: oPUMAS.jl after M5 (2026-10-02)

- 🔖 **Commit:** `a2abde0` ✨ feat: 🎛️⚖️ M5 switch matrix and physical checks of oPUMAS.jl (branch `dev/julia`, pushed to `origin`).
- ✅ **What was done:**
  - 🔁 `julia/test/test_step.jl` runs the full step for all 16 switch rows that have a golden file. Every row is bit-for-bit: 155 of 155 fields, and 159 of 159 for `sb2001`.
  - ⚖️ The new `julia/test/test_physics.jl` adds four checks:
    - water budget per column;
    - energy budget per column, with ice as the reference phase;
    - negative mixing ratios, with the allowed margin set by PUMAS' `qsmall Δt` limiter;
    - a clear dry column must have zero tendencies.
  - ✅ `Pkg.test()`: 2751 of 2751 pass, about 90 s on herdr-compute.
- 📌 **Current state:**
  - 🏁 M0–M5 are done.
  - ⏭️ `tau` and `wsub_default` have no golden file and are skipped.
  - 🚫 The `error(...)` stubs (Bigg freezing, evaporation without graupel) are never reached by the in-scope rows.
- ⏭️ **Next steps:**
  1. 🌿 Merge `dev/julia` into `main_cam_interface`, per the plan. Waiting for Diego's go-ahead. Then run `translation/sync.sh` and `python3 translation/check.py --branches`.
  2. 📖 M6 readability pass: replace the PUMAS short names in `src/stages/` with container names, keeping the tests bit-for-bit.
  3. 🔍 Optional:
     - trace the NaN in CAM `tnd_nsnow`/`re_ice`;
     - time the Fortran driver;
     - cut the 17.5 MiB of allocations per step.
- 🧭 **Key decisions:**
  - ➖ Negatives are allowed down to `qsmall Δt` plus 1e-12 of column activity. The rain and snow limiters only act on sink rates ≥ `qsmall` (FINDINGS.md), so this bound comes from PUMAS, not from the twin.
  - 📏 Budget tolerances are 1e-12 × gross terms, plus a floor of 1e-17 kg m⁻² s⁻¹ for water and 1e-11 W m⁻² for energy.
