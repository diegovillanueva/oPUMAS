# oPUMAS

🧊🌧️ Diego's private development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) (Parameterization for Unified Microphysics Across Scales).
🔄🌿 It tracks upstream fixes continuously; changes that should go upstream will later be sent as a single PR to ESCOMP/PUMAS.

- 🧭🧊 What separates oPUMAS from PUMAS (naming, defaults, sync): [`PHILOSOPHY.md`](PHILOSOPHY.md); name dictionary, translation scripts and the `CF-CCPP` branch: [`translation/README.md`](translation/README.md)
- 📘📄 Original upstream README (checkout, running with CAM): [`README_PUMAS.md`](README_PUMAS.md)
- 🐙🔁 Git setup, daily upstream sync, and (disabled) PR recipe: [`git_work.md`](git_work.md)
- 🤖📋 Agent rules: [`AGENTS.md`](AGENTS.md)
- 🧪♻️ Compile check, bit-for-bit test and switch matrix: [`tests/README.md`](tests/README.md)
- 🤝📂 Handoffs of all topics: [`handoffs/`](handoffs/)
- 🧊🔌 Run oPUMAS inside CESM 2.1.5 CAM (`micro_mg_version = 3`): [`interface/cesm215/README.md`](interface/cesm215/README.md)

## Remotes and branches

| name | what | rule |
|------|------|------|
| `origin` | github.com/diegovillanueva/oPUMAS (private) | push freely |
| `upstream` | github.com/ESCOMP/PUMAS | fetch only (push URL = `DISABLED`) |
| `main_cam` | mirror of `upstream/main_cam` | fast-forward only, never commit |
| `main_cam_interface` | Diego's development in PUMAS names (default branch) | merge `main_cam` in regularly; dictionary edits happen here |
| `CF-CCPP` | the same code in CF-CCPP names, editable | run `translation/sync.sh` before and after working on it |

## Local clone

```
/cluster/work/climate/dvillanueva/work/CESM/PUMAS
```

🌿📂 Worktree of `CF-CCPP` (complete `-usr_src` directory for CESM 2.1.5 under `interface/cesm215`):

```
/cluster/work/climate/dvillanueva/work/CESM/cesm215_pumas_worktrees/opumas_CF-CCPP
```

## ✅🧪 Validate

```bash
tests/compile_check.sh                    # PUMAS files and CAM overrides compile (about 40 s, herdr-compute)
tests/bfb.sh run                          # 1-day CESM 2.1.5 run is bit-for-bit with the baseline (about 15 min)
python3 -m pytest -q translation          # translator tests
python3 translation/check.py              # dictionary against the sources
python3 translation/check.py --branches   # CF-CCPP is exactly the translation of main_cam_interface
```

🧊♻️ Bit-for-bit test in CESM 2.1.5 (cases `bfb_dev` and `bfb_cfccpp`): see [`interface/cesm215/README.md`](interface/cesm215/README.md).

## 🔬📂 Code analysis notes

- 📝 `code_analysis/` holds per-topic notes on where schemes live in the source (e.g. `code_analysis/PUMAS_immersion_freezing.md`).
- 🗺️ Start with the stage map of `micro_pumas_tend`: [`code_analysis/PUMAS_stage_map.md`](code_analysis/PUMAS_stage_map.md) (16 stages with line ranges, search anchors and a checklist for adding a process).
