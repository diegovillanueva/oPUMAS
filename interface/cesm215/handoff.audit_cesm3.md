# 🤝🔎 Handoff: audit of the cesm215 interface against CESM 2.3 / CESM 3 (2026-10-01)

- 📌 Commit: see `git log -1 -- interface/cesm215/README.md` (docs-only commit on `dev`).
- ✅ What was done: compared graupel handling of the interface with `cam6_3_097` (local) and `cam6_4_204` (GitHub, CAM7/CESM3-beta physics). Result: graupel is prognostic by default in CAM7 but `graupel_in_rad` defaults to `.false.` everywhere, so rimed ice is **not** radiatively active in CESM 3 beta either. The interface already matches; no code change needed.
- 🧊 Current state: interface = rel_60 host + PUMAS v1.41 + switchable graupel (off by default), graupel invisible to radiation, 1-day graupel run validated.
- ➡️ Next steps (optional): (1) add `graupel_in_rad`-style switch with `DEGRAU`/`ICGRAUWP`/`CLDFGRAU` pbuf fields and overrides of `cloud_rad_props.F90`/`radiation.F90` ported from cam6_3_097; (2) add the missing graupel diagnostics (`LS_REFFGRAU`, `CLDFGRAU`, `UMG`, `FREQG`, `QGSEDTEN`, graupel process rates); (3) 3-month `pumas_grau` run.
- 🧭 Key decisions: keep cam6 defaults rather than CAM7 MG3 tuning (`vtrmi 1.39`, `accre_enhan 0.731`, `implicit_fall`, `accre_sees_auto`, `precip_fall_corr`), documented in the README so they can be set via `user_nl_cam`.
