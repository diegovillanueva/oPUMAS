# 🤝🌍 Handoff: CESM 3-like namelist and activation-velocity switches (2026-10-01)

- 📌 Commit: see `git log -1 -- interface/cesm215/user_nl_cam.pumas_cesm3like`.
- ✅ What was done: new overrides `microp_aero.F90` and `ndrop.F90` (CAM7 `microp_aero_wsub_scale`, `wsubi_scale`, `wsub_min`, `wsub_min_asf`, `wsubi_min`; defaults = rel_60 hard-coded values); XML entries committed on the CAM worktree branch `dv_pumas_rel60` (09f11e59) and in the regenerated patch; `user_nl_cam.pumas_cesm3like`; `TODO.md`.
- 📌🔖 Reference for "CESM 3": ESCOMP/CESM `cesm3_0_beta10` -> ESCOMP/CAM `cam6_4_203` (micro defaults identical in `cam6_4_204`), `atmos_phys0_29_000` for `scale_subgrid_vertical_velocity`.
- 🧪 Current state: (1) default namelist is bit-for-bit: `pumas_grau` rebuilt with the new overrides gives a byte-identical `cam.r` after 1 day. (2) `cesm215_pumas.F2000climo.pumas_c3like` (clone of `pumas_grau`, `--keepexe`) ran 1 day with the CESM 3-like namelist; fields differ from `pumas_grau` as expected (global unweighted means: GRAUQM 1.3e-8 vs 4.6e-9, CLDLIQ 3.3e-6 vs 3.9e-6, NUMLIQ 9.5e5 vs 1.1e6).
- ⚠️ Gotcha: `create_clone` resets `SRCROOT` to `cesm215/cime/..`; fix it in `env_case.xml` and `LockedFiles/env_case.xml` to `$W/cesm215_pumas_worktrees/cesm215_pumas`, otherwise build-namelist rejects the oPUMAS variables.
- ➡️ Next steps: see `TODO.md`.
- 🧭 Key decisions: `wsub_min_asf` is implemented as a public module variable `wmixmin` in `ndrop.F90` (set by `microp_aero_readnl`) instead of a new `dropmixnuc` argument, so SPCAM callers stay untouched. `microp_aero_npccn_scale` was not ported (CAM7 default 1).
