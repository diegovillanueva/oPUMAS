# 💡🧩 Proposal: replace the full-file CAM overrides

📌🚧 Status: proposal only, nothing here is implemented.

## 🧱⚠️ Problem

- 📄🔁 Five CAM files are stored as full copies (7645 lines) to carry 465 changed lines; `micro_mg_cam.F90` alone accounts for 379 of them.
- 🔒🏷️ The copies are frozen at `cam_cesm2_1_rel_60`. Porting to another tag or to `dv_diag` means re-deriving the patch by hand.
- 🔎🙈 A review of the repo shows thousands of lines of CAM code that are not ours; the real change is only visible through `make_interface.sh diff`.

## ✅🎯 Recommended design: one new glue module plus a small patch series

1. 🧊🔌 **Move everything PUMAS-specific into a new file `micro_pumas_cam.F90`** in this directory, as upstream CAM does. It holds the `micro_pumas_nl` namelist read, the `micro_pumas_init` call, the packing and the `micro_pumas_tend` call. It is a new module, not an override, so `-usr_src` simply adds it. It is ours from the first line.
2. 🩹📏 **Reduce the CAM edits to hooks.** `micro_mg_cam.F90` keeps only a `case (3)` that calls the new module in readnl, init, get_cols and tend. Together with the small edits in `wv_sat_methods`, `check_energy`, `microp_aero` and `ndrop` this is on the order of 100 lines that contain no PUMAS names.
3. 📦🌿 **Store those hooks as a patch series** in `patches/000N-*.patch` and apply them with `git am` to the CAM worktree branch `dv_pumas_rel60`, exactly as the namelist XML patch is handled today. The five full copies leave the repo, and `-usr_src` contains only PUMAS files and the glue module.
4. 🛠️🔁 **Extend `make_interface.sh`** with `export` (regenerate the patches from the worktree branch) and `verify` (fail if the branch and the patches differ).

## ⚖️📋 What it buys and what it costs

| | Today | Proposed |
|---|---|---|
| Lines of foreign CAM code in the repo | about 7200 | 0 |
| Where the real change is visible | `make_interface.sh diff` | the repo itself |
| Port to a new CAM tag or `dv_diag` | manual | `git am -3` or rebase of a 100-line series |
| Editing a CAM hook | edit the copy | edit the worktree, then `make_interface.sh export` |
| Zero-edit promise for the CAM tree | XML only | XML plus hooks, on a dedicated branch |

- ⚠️🧶 **Main risk.** The `case (3)` code uses module-private variables of `micro_mg_cam` (pbuf indices, constituent indices, history flags). The glue module needs them passed as arguments or made public, and this decides whether step 1 is clean.
- 🪶🔀 **Cheaper alternative.** Skip step 1 and only do steps 3 and 4. The repo then holds a 465-line patch instead of 7645 lines, but the PUMAS glue stays inside a CAM file.

## 🪜🧪 Migration, each step verified bit-for-bit

1. 📦 Steps 3 and 4 first: same compiled source, so `tests/bfb.sh run` must pass without a new baseline.
2. 🧊 Then step 1 in small moves (readnl, init, tend), each with `tests/compile_check.sh` and `tests/bfb.sh run`.
3. 🎚️ Before step 1, fill the baselines of `tests/switches.csv`, so every switch path is protected.
