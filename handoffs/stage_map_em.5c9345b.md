# 🤝🗺️ Handoff: Vim-friendly stage map

- 🔖 **Commit**: `5c9345b` 📝 docs: 🖥️👓 Vim-friendly emoji stage map (PUMAS_stage_map.em.md)
- ✅📝 **What was done**: added `code_analysis/PUMAS_stage_map.em.md`, a table-free, emoji-coded twin of `PUMAS_stage_map.md`, linked from `README.md`.
- 📍🧊 **Current state**: both versions carry line numbers of `main_cam_interface:micro_pumas_v1.F90`; on `dev/Sppe` they are shifted by +4..+37 lines.
- ⏭️🔁 **Next steps**: when the table version changes, mirror the edit in the `.em.md` file (or refresh line numbers after an upstream merge).
- 🧠🎨 **Key decisions**: one block per stage instead of a table, colour squares for stage kind (🟦 setup, 🟩 rates, 🟥 direct state change, 🟨 tendencies, 🟪 diagnostics), lines kept ≤ 100 chars.
