# oPUMAS

🧊🌧️ Diego's private development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) (Parameterization for Unified Microphysics Across Scales).
🔄🌿 It tracks upstream fixes continuously; changes that should go upstream will later be sent as a single PR to ESCOMP/PUMAS.

- 🧭🧊 What separates oPUMAS from PUMAS (naming, defaults, sync): [`PHILOSOPHY.md`](PHILOSOPHY.md); name dictionary: [`name_map.csv`](name_map.csv)
- 📘📄 Original upstream README (checkout, running with CAM): [`README_PUMAS.md`](README_PUMAS.md)
- 🐙🔁 Git setup, daily upstream sync, and (disabled) PR recipe: [`git_work.md`](git_work.md)
- 🤖📋 Agent rules: [`AGENTS.md`](AGENTS.md)
- 🧊🔌 Run oPUMAS inside CESM 2.1.5 CAM (`micro_mg_version = 3`): [`interface/cesm215/README.md`](interface/cesm215/README.md)

## Remotes and branches

| name | what | rule |
|------|------|------|
| `origin` | github.com/diegovillanueva/oPUMAS (private) | push freely |
| `upstream` | github.com/ESCOMP/PUMAS | fetch only (push URL = `DISABLED`) |
| `main_cam` | mirror of `upstream/main_cam` | fast-forward only, never commit |
| `dev` | Diego's development (default branch) | merge `main_cam` in regularly |

## Local clone

```
/cluster/work/climate/dvillanueva/work/CESM/PUMAS
```

## 🔬📂 Code analysis notes

- 📝 `code_analysis/` holds per-topic notes on where schemes live in the source (e.g. `code_analysis/PUMAS_immersion_freezing.md`).
