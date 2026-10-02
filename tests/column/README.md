# 🧪📦 Column golden data for the Julia twin

🎯 Reference inputs and outputs of `micro_pumas_tend`, taken from CESM 2.1.5 and from a standalone Fortran driver, so `julia/` can be tested column by column and process by process.

## 🗺️🧭 Pieces

| File | What it does |
|---|---|
| `driver.F90` | 🛠️ Standalone program: reads `settings.nml` and an inputs file, calls `micro_pumas_init` and `micro_pumas_tend` exactly as `interface/cesm215/micro_mg_cam.F90` does, writes all outputs and every allocated `proc_rates` field. |
| `pumas_column_io.F90` | 📤 Shared NetCDF writer used by the driver and by the CAM dump, so both files have the same layout. |
| `build.sh` | 🧱 Builds the driver with gfortran 12.2 at `-O0 -ffp-contract=off` (about 30 s). `OPT=-O2` mimics the CESM build; `CAPTURE=1` builds the capture variant. |
| `capture/` | 🔬 `patch_v1.py` copies `micro_pumas_v1.F90` with a record of every process call (76 call sites) before and after the call; `pumas_column_capture.F90` writes them. |
| `dump/` | 📥 `make_case.sh` clones `bfb_dev` into `cesm215_pumas.F2000climo.column_dump` with a `-usr_src` copy of `interface/cesm215` whose `micro_mg_cam.F90` is patched by `patch_micro_mg_cam.py`; `run_case.sh` runs one day and dumps chosen steps and chunks. |
| `select_columns.py` | 🔎 Picks about 200 dumped columns so every process rate is non-zero in at least 5 of them, plus 15 per cloud regime. |
| `synthetic.py` | 🌦️ Nine hand-made columns for rare processes (supercooled rain, homogeneous freezing, melting, dry evaporation, riming, drizzle, ice supersaturation, warm ice, clear sky). |
| `compare.py` | ⚖️ Variable-by-variable comparison with `--rtol`/`--atol`; zero tolerances mean bit-for-bit. |
| `run.sh` | ▶️ Rebuilds everything below from the stored columns in under a minute. |

## ▶️🧾 How to run

```bash
# 1️⃣ once per CESM dump (only after a change in micro_mg_cam.F90 or the case):
tests/column/dump/make_case.sh          # clone, patch, build (sbatch --wait, about 25 min)
tests/column/dump/run_case.sh           # 1-day run, dumps steps 24 and 47, every 50th chunk
tests/column/run.sh select /cluster/work/climate/dvillanueva/cesm_runs/cesm215_pumas.F2000climo.column_dump/run/column_dump
# 2️⃣ every time PUMAS changes (inside herdr-compute):
tests/column/run.sh
```

📂 Stored columns (`SOURCE`): `/cluster/work/climate/dvillanueva/opumas_column/source/` with `inputs.nc`, `cam_outputs.nc`, `synthetic_inputs.nc`, `settings.nml` and `coverage.csv`.
📂 Golden files (`DATA`, git-ignored): `tests/column/data/` with `golden_inputs.nc`, `golden_<label>.nc` per switch row, `golden_processes.nc` (base settings), `golden_processes_sb2001.nc` (Seifert-Beheng warm rain, whose routines never run in the base settings) and `selfcheck.log`.

## ✅🔁 Checks done by `run.sh`

- ♻️ **Self-check.** The driver, fed CAM's own inputs and settings, reproduces CAM's outputs for the dumped columns (`selfcheck.log`). This proves the driver calls PUMAS the way CAM does.
- 🔬 **Capture is neutral.** The capture build gives the same outputs as the plain build, bit for bit.

## 📊🔎 Coverage of the stored columns

📥 Dump of 2026-10-02: case `cesm215_pumas.F2000climo.column_dump` (bfb_dev settings), steps 24 and 47, every 50th chunk, 2139 columns. `select_columns.py` keeps 200 of them; `coverage.csv` lists every rate.

| Rate | Selected columns with rate ≠ 0 | Dumped columns with rate ≠ 0 |
|---|---|---|
| all 7 cloud regimes | 16 to 179 | 53 to 1824 |
| deposition nucleation `mnudeptot` | 5 | 12 |
| contact freezing `mnuccttot` | 5 | 20 |
| rime splintering on rain-graupel `qmultrgtot` | 6 | 15 |
| rime splintering on graupel `qmultgtot` | 7 | 15 |
| heterogeneous rain freezing `mnuccrtot` | 9 | 24 |
| `qcrestot`, `qirestot`, `qvres`, `qcsevap`, `qisevap` | 0 | 0 |

🚫 The last row is zero by construction. `qcrestot`, `qirestot` and `qvres` need `remove_supersat`, which CAM fills from `allow_sed_supersat = .false.`. `qcsevap` and `qisevap` are only computed by the explicit `Sedimentation`, and these settings use implicit fall; the `explicit_fall` switch row covers them.
✅ Self-check on this dump: 155 variables pass, worst relative error 2.4e-16 (`-O0`, `SELF_RTOL=1e-12`).

## 🎛️🚫 Switch rows

🎛️ Every row of `tests/switches.csv` gets `golden_<label>.nc`, made from `settings.nml` plus the row's namelist overrides.
🧊🎛️ Sppe rows: `switches_sppe.csv` lists the 11 process controls of `dev/Sppe` (`micro_pumas_sppe.F90`). `run.sh` builds a second driver with `SPPE=1 build.sh` (into `build/sppe`), which takes `micro_pumas_sppe.F90`, `micro_pumas_utils.F90` and `micro_pumas_v1.F90` from `git archive dev/Sppe` (`SPPE_REF` to change it) without merging that branch, and checks that `golden_sppe_base.nc` is bit-for-bit with `golden_base.nc`. `SPPE_ROWS=0 run.sh` skips them. `sppe_no_cnt_frz` equals base by design: `onoff_contact_freezing` only acts on the Bigg path used without host freezing, which the Julia twin does not translate.

🌧️🧮 The `tau` row reads `KBARF_tau_kernel.dat` from the repo root (`@REPO@` in `switches.csv` becomes the repo path in `settings_tau.nml`).

🚫 Skipped row: `wsub_default` (changes `microp_aero`, not PUMAS, so it has no meaning for fixed inputs).

## 📐🧷 File layout notes

- 🔄 Files are written by Fortran, so xarray and NCDatasets.jl see `(lev, col)` in Python order and `(col, lev)` in Julia order. Interface levels use `ilev`, dust bins `dust_bin`.
- 🔬 In `golden_processes.nc` each variable is `<routine>_<n>__<argument>__<in|out>`. `<n>` counts call sites of the routine in source order, `<argument>` is `aNN` for the NN-th positional argument or the keyword name, and the attribute `expr` holds the Fortran actual argument. A second execution of the same site gets `_v2`. 1-D vectors over all packed points are stored as `(col, lev)`.
