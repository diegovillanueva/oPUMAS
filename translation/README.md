# 🔤📖 Name dictionary: PUMAS Fortran ↔ oPUMAS.jl

🎯📘 The Fortran keeps the upstream PUMAS names, so upstream merges stay easy. The long, readable names live in the Julia twin [`julia/`](../julia/README.md). [`dictionary.csv`](dictionary.csv) maps one to the other, one row per concept.
🗄️🌿 Until 2026-10-03 a translator script kept a second Fortran branch `CF-CCPP` in long names. It was retired with the move to a single `main` branch; its last state is the tag `archive/CF-CCPP`, and the script is in the git history (`git show archive/CF-CCPP:translation/translate.py`).

## 📖🧱 Columns

| Column | Content |
|---|---|
| `id` | 🔑 Unique key: the Fortran name, `proc_rates%<name>` for a component, `namelist:<name>` for a setting |
| `kind` | `rate`, `state`, `size`, `fallspeed`, `sedimentation`, `air`, `fraction`, `saturation`, `tendency`, `diagnostic`, `flux`, `setting` |
| `scope` | `global` (local variable of the core), `type:proc_rates_type` (component), `argument:micro_pumas_tend` (output argument), `namelist` |
| `units` | 🧮 Udunits-like (`kg kg-1 s-1`); `1` or empty for dimensionless |
| `meaning` | 📝 One line for humans |
| `Fortran` | 🧊 The name as spelled in upstream PUMAS (for a setting: the namelist variable) |
| `Julia` | 🟣 The path in oPUMAS.jl: a workspace field (`rates.rain.accretion_by_snow.mass`, `diagnostics.surface.rain_rate`) or a `Settings` field (`settings.splinter_number_factor`) |
| `standard_name` | 📚 The CF/CCPP standard name when one exists |
| `notes` | 🗒️ Free text |

## 🔁✅ Keeping dictionary, Julia and Fortran in sync

🧪🔗 `julia/test/names_check.jl` runs inside `Pkg.test()` (test set M7) and fails when any of these break:

1. 🧾 The header is the one above, ids are unique, and no `Fortran`, `Julia` or `meaning` cell is empty.
2. 🟣 Every `Julia` path exists in the workspace built by `allocate`, or as a field of `Settings`.
3. 🔁 Every entry of the name tables in `julia/src/names.jl` and of `NAMELIST_FIELDS` in `julia/src/namelist.jl` has its row (the other direction).
4. 🧊 Every `Fortran` name still occurs in the repo's Fortran sources, which catches upstream renames after a merge of `main_cam`.
5. 📜 Every exported Julia function names its Fortran routine in its docstring, as listed in `julia/TWINS.csv`.

```bash
julia/tools/julia.sh -e 'using Pkg; Pkg.test()'                        # all Julia tests, M7 included
julia/tools/julia.sh --project=julia julia/tools/check_names.jl        # only the name checks (about 10 s)
```

## ➕📖 Add or change a name

1. 🟣 Add the field in oPUMAS.jl (and its entry in `src/names.jl` or `NAMELIST_FIELDS` when it maps a Fortran name).
2. 📖 Add or edit the row in `dictionary.csv` in the same commit.
3. ✅ Run `check_names.jl` or `Pkg.test()`; fix every error.

## ✍️📐 Naming grammar of the Julia paths

```
rates          rates.<species>.<process>[_<by|to|from>_<partner>].<mass|number>   local (in-cloud) rates
process_rates  process_rates.<species>.<process>...                              grid means (proc_rates)
outputs        tendencies.<...>, diagnostics.<...>, diagnostics.flux.<...>, diagnostics.surface.<...>
settings       settings.<field> of Settings (namelist values)
```

- 🔎📚 Words come from CCPP and CF standard names, then the AMS Glossary (autoconversion, accretion, riming, aggregation, deposition, sublimation).
- 🏷️🧊 Species: `cloud_liquid`, `cloud_ice`, `rain`, `snow`, `graupel`, `vapor`.
- 🧮🚫 No unit suffix in Julia names; units are in the `units` column and in docstrings.
- 🔒 Names that hosts see (module names, namelist variables, history fields, `micro_pumas_init`/`micro_pumas_tend` and their arguments) are never renamed in the Fortran.
- 🧮 Not yet in the dictionary: Fortran scratch variables (`dum*`, `t`, `q`), TAU bin arrays and fixer diagnostics.
