# 🔤🌍 translation/: name dictionary and branch sync

🎯📘 oPUMAS keeps the PUMAS identifiers on `dev` and offers the same code with self-explaining names on the branch `CF-CCPP`. Both branches are editable. One dictionary and one script keep them equivalent.

| File | Role |
|---|---|
| `dictionary.csv` | 📖 One row per concept, one column per language (`PUMAS`, `CF-CCPP`, ...) |
| `languages.csv` | 🌿 Language, its branch, and whether unit suffixes are appended |
| `files.csv` | 📂 Files that are translated and how (`full` or `components`) |
| `abbreviations.csv` | ✂️ The only shortenings allowed in constructed names |
| `translate.py` | 🔁 Translates identifiers between any two languages |
| `check.py` | ✅ Static checks of the dictionary and the branch invariant |
| `sync.sh` | 🔀 Two-way sync between `dev` and the language branches |
| `test_translate.py` | 🧪 Pytest on small Fortran snippets and on the real dictionary |

## 🧭🔧 Daily use

```bash
translation/sync.sh                      # before and after working on CF-CCPP
python3 translation/check.py             # dictionary against the sources
python3 translation/check.py --branches  # translate(dev) == CF-CCPP
python3 -m pytest -q translation
git push origin dev CF-CCPP              # sync.sh never pushes
```

- 🌿🔁 `sync.sh` replays new `CF-CCPP` commits, translated, onto `dev`, then writes one `🔁 sync` commit on `CF-CCPP` that holds the translated `dev` tip. The `Translated-from: <hash>` trailer of that commit marks the synced pair.
- 📥🧬 Upstream PUMAS is merged into `dev` only. Dictionary edits are made on `dev` only. The next sync carries both to `CF-CCPP`.
- 🧳📤 One-way export for other compilers: `translate.py --from PUMAS --to CF-CCPP --wrap` breaks lines over 132 characters, and `--no-units` drops the unit suffixes. Never use these two flags on the editable branches.

## 📖🧱 Dictionary columns

| Column | Content |
|---|---|
| `id` | 🔑 Unique key: the PUMAS name, or `proc_rates%<name>` for a component |
| `kind` | `rate`, `state`, `size`, `fallspeed`, `sedimentation`, `air`, `fraction`, `saturation`, `diagnostic` |
| `scope` | `global` (all `full` files), `type:proc_rates_type` (component), or `file:<name.F90>` |
| `units` | 🧮 Udunits-like (`kg kg-1 s-1`); the suffix (`_kg_per_kg_per_s`) is generated from it, `1` or empty gives none |
| `meaning` | 📝 One line for humans |
| `PUMAS`, `CF-CCPP` | 🔤 The name in each language, without unit suffix |
| `cf_ccpp_source` | `CCPP`, `CF` or `constructed` |
| `standard_name` | 📚 The strict CF/CCPP standard name when one exists |
| `notes` | 🗒️ Free text |

## ✍️📐 Naming grammar of the CF-CCPP language

```
rates        <species>_<mass|number>_tend_<process>[_<by|to|from>_<partner>]
state        <species>_mass_mixing_ratio, <species>_number_conc, incloud_/inprecip_ prefix
sizes        <species>_size_dist_<slope|intercept|shape>
fall speeds  <species>_<mass|number>_weighted_fall_speed
```

- 🔎📚 Lookup order: CCPP standard name, then CF standard name, then a constructed name with AMS Glossary process words (autoconversion, accretion, riming, aggregation, deposition, sublimation).
- 🏷️🧊 Species words: `cldliq`, `cldice`, `rain`, `snow`, `graupel`, `vapor`.
- 🌫️📦 `proc_rates` components are grid means and carry the plain name. Their local twins in `micro_pumas_tend` are in-cloud or in-precipitation rates and carry the prefix `local_`.
- 📏⛔ A name plus its unit suffix has at most 63 characters (Fortran limit). `check.py` enforces it.

## ➕📖 Add or change a name

1. 🌿 On `dev`, add or edit a row of `dictionary.csv`.
2. ✅ Run `python3 translation/check.py`; fix every error.
3. 💾 Commit on `dev`, run `translation/sync.sh`, push both branches.

## ➕🌐 Add a language

1. 📖 Add a column to `dictionary.csv` (any header that is not a metadata column). Empty cells mean "keep the PUMAS name".
2. 🌿 Add a row to `languages.csv` with the branch name.
3. 🔁 Run `translation/sync.sh`; it creates the branch.

## 🛡️⚠️ Rules the checker enforces, and known limits

- 🔒 Public names are never translated: module names, `micro_pumas_init`, `micro_pumas_tend`, their dummy arguments, names imported from untranslated modules, namelist and history names.
- 🔤 Fortran is case-insensitive. A PUMAS name must be spelled as in the dictionary or all upper case. `dv` is not in the dictionary because the source also spells it `Dv`.
- 💬 Comments, strings and preprocessor lines are never touched; `!$acc` and `!$omp` lines are translated.
- 🌍 A `global` row renames the identifier in every `full` file, also where a routine uses the same short name for a slightly different quantity (for example `qc` inside `pumas_stochastic_collect_tau.F90`). This is always safe for the compiler; the name can be less exact there.
- 📏 Translated lines exceed 132 characters. The CESM gnu build uses unlimited line length, as upstream PUMAS already needs it.
- 🧮 Not yet in the dictionary: scratch variables (`dum*`, `t`, `q`), TAU bin arrays, fixer and scaling diagnostics, about 600 further identifiers.
