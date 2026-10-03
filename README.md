# oPUMAS

🧊🌧️ Diego's private development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) (Parameterization for Unified Microphysics Across Scales).
🔄🌿 It tracks upstream fixes continuously; changes that should go upstream will later be sent as a single PR to ESCOMP/PUMAS.

- 🧭🧊 What separates oPUMAS from PUMAS (naming, defaults, sync): [`PHILOSOPHY.md`](PHILOSOPHY.md); Fortran ↔ Julia name dictionary: [`translation/README.md`](translation/README.md)
- 📘📄 Original upstream README (checkout, running with CAM): [`README_PUMAS.md`](README_PUMAS.md)
- 🐙🔁 Git setup, daily upstream sync, and (disabled) PR recipe: [`git_work.md`](git_work.md)
- 🤖📋 Agent rules: [`AGENTS.md`](AGENTS.md)
- 🧪♻️ Compile check, bit-for-bit test and switch matrix: [`tests/README.md`](tests/README.md); column driver and golden data: [`tests/column/README.md`](tests/column/README.md)
- 🤝📂 Handoffs of all topics: [`handoffs/`](handoffs/)
- 🟣📖 Readable Julia twin (in progress): [`julia/README.md`](julia/README.md); plan in [`handoffs/julia.plan.md`](handoffs/julia.plan.md)
- 🧊🔌 Run oPUMAS inside CESM 2.1.5 CAM (`micro_mg_version = 3`): [`interface/cesm215/README.md`](interface/cesm215/README.md)
- 🌍🩹 Sppe CAM-side patches for CESM 3 / CAM7 (checked, untested): [`interface/cesm3/README.md`](interface/cesm3/README.md)

## Remotes and branches

| name | what | rule |
|------|------|------|
| `origin` | github.com/diegovillanueva/oPUMAS (private) | push freely |
| `upstream` | github.com/ESCOMP/PUMAS | fetch only (push URL = `DISABLED`) |
| `main_cam` | mirror of `upstream/main_cam` | fast-forward only, never commit |
| `main` | Diego's development (default branch): upstream PUMAS plus Sppe controls (`micro_pumas_sppe.F90`, group `micro_dv_nl`, all off by default), tests, dictionary and the Julia twin | merge `main_cam` in regularly; `feat/*` branches merge back here |
| `archive/*` (tags) | `archive/CF-CCPP` (retired long-name Fortran branch), `archive/dev-Sppe`, `archive/dev-julia` | read only |

🧪🏷️ `bfb_dev` in `tests/` is a CESM case name, not a branch. 🗓️🌿 Until 2026-10-03 the development branch was called `main_cam_interface`.

## Local clone

```
/cluster/work/climate/dvillanueva/work/CESM/PUMAS
```

## ✅🧪 Validate

```bash
tests/compile_check.sh                    # PUMAS files and CAM overrides compile (about 40 s, herdr-compute)
tests/bfb.sh run                          # 3-step CESM 2.1.5 run is bit-for-bit with the baseline
tests/bfb.sh par all                      # every switch row in parallel, one --keepexe clone each
tests/column/run.sh                       # column driver self-check, golden files and 3-day loop for julia/ (about 1 min)
julia/tools/julia.sh -e 'using Pkg; Pkg.test()'   # Julia twin against the golden files, dictionary in sync
```

🧊♻️ Bit-for-bit test in CESM 2.1.5 (case `bfb_dev`): see [`interface/cesm215/README.md`](interface/cesm215/README.md).

## 🔬📂 Code analysis notes

- 📝 `code_analysis/` holds per-topic notes on where schemes live in the source (e.g. `code_analysis/PUMAS_immersion_freezing.md`).
- 🗺️ Start with the stage map of `micro_pumas_tend`: [`code_analysis/PUMAS_stage_map.md`](code_analysis/PUMAS_stage_map.md) (16 stages with line ranges, search anchors and a checklist for adding a process).
- 🖥️👓 Vim-friendly emoji twin without tables: [`code_analysis/PUMAS_stage_map.em.md`](code_analysis/PUMAS_stage_map.em.md) (`*.em.md` = plain-text layout, readable in a terminal).
