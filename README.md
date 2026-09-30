# oPUMAS

🧊🌧️ Diego's private development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) (Parameterization for Unified Microphysics Across Scales).
🔄🌿 It tracks upstream fixes continuously; changes that should go upstream will later be sent as a single PR to ESCOMP/PUMAS.

- 📘📄 Original upstream README (checkout, running with CAM): [`README_PUMAS.md`](README_PUMAS.md)
- 🐙🔁 Git setup, daily upstream sync, and (disabled) PR recipe: [`git_work.md`](git_work.md)
- 🤖📋 Agent rules: [`AGENTS.md`](AGENTS.md)

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
