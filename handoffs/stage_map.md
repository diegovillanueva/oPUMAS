# 🤝🗺️ Handoff: stage map of micro_pumas_tend, deferred items in TODO.md (2026-10-01)

- 📌 Commit: see `git log -1 -- code_analysis/PUMAS_stage_map.md`.
- ✅ What was done: `code_analysis/PUMAS_stage_map.md` (16 stages, line ranges, search anchors, checklist for adding a process). `interface/cesm215/TODO.md` gained two deferred items: fill the switch baselines (with the proposed test loop) and replace the full-file overrides (four-step plan).
- 🧭 Key decisions: Diego chose the stage map instead of a structural refactor of the core, so `dev` stays identical to upstream. Issue 1 keeps the structural split only as an option. The three-month runs and missing interface features are out of scope for now.
- 🧪 Current state: all 16 search anchors were checked to exist in `micro_pumas_v1.F90` at the stated lines. The stage descriptions come from reading comments and calls, not from running the code. Docs only, no Fortran changed.
- ➡️ Next steps: the two deferred items in `TODO.md`, baselines first.
