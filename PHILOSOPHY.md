# 🧭🧊 oPUMAS philosophy: what separates oPUMAS from PUMAS

🌱📘 oPUMAS is a development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) that aims to be easier to read, learn and extend, without drifting away from PUMAS.
🧩📐 Every design decision in this repo is checked against the principles below; a new principle is added here before it is applied in the code.

## 1. 🔤📖 Names are for readers, the dictionary is for the sync

- 🔤✨ oPUMAS offers variable and function names that minimise the effort of understanding what they mean.
- 📖🔁 Every such name keeps a dictionary entry that maps it to its reference name in PUMAS, so that syncing PUMAS commits into oPUMAS (and sending changes back) stays easy.
- 🌿🌿 The same code lives on two editable branches: `main_cam_interface` in the PUMAS language and `CF-CCPP` in a language built from CCPP and CF standard names and AMS Glossary process words. `translation/sync.sh` keeps them equivalent, so everybody works in the language they read best.
- 📄🗂️ The dictionary is [`translation/dictionary.csv`](translation/dictionary.csv): one row per concept, one column per language. A new language is a new column plus a row in `translation/languages.csv`. Details: [`translation/README.md`](translation/README.md).
- 🧮🏷️ Units are part of the name: the dictionary has a `units` column and the translator appends the suffix (`npracs` becomes `local_rain_number_tend_accretion_by_snow_per_kg_per_s`). `--no-units` switches the suffix off.
- 🤖🔄 The translation is done by a script (`translation/translate.py`), never by hand, in both directions.
- ⚖️🧪 A translation never changes results: `CF-CCPP` is verified bit-for-bit against `main_cam_interface` in CESM 2.1.5, and `translation/check.py --branches` verifies that `CF-CCPP` is exactly the translation of `main_cam_interface`.
- 🎛️🤝 Names that users and host models see (module names, namelist variables, history field names, `micro_pumas_init`/`micro_pumas_tend` and their argument lists) keep their PUMAS spelling in every language, so existing namelists, scripts and host interfaces keep working. The components of `proc_rates_type` are translated, together with the host files that read them.
- 🚧📅 Status (2026-10-01): 219 names are in the dictionary (process rates, hydrometeor state, size-distribution, fall-speed and sedimentation names). Scratch variables and about 600 further identifiers still carry their PUMAS names on both branches.

## 2. 🎛️🔕 Everything new is switchable and defaults to PUMAS

- 🎚️✅ Every feature added in oPUMAS or in an `interface/*` is switchable and off by default: an empty namelist makes oPUMAS behave exactly like upstream PUMAS, bit-for-bit where possible.
- 🧭📌 When out-of-the-box CESM–PUMAS versions differ, the CESM 2.3 behaviour is the default; settings of other versions are offered as namelist files (e.g. `interface/cesm215/user_nl_cam.pumas_cesm3like`).

## 3. 🔄🌿 Stay close to upstream

- 🔁📥 Upstream PUMAS commits are merged continuously (`main_cam` mirrors `upstream/main_cam`; see [`git_work.md`](git_work.md)).
- 🧱🔌 Host-model glue lives in `interface/<host>/`, never in the oPUMAS core files, so the core stays comparable with PUMAS.

## 4. 📚🪜 Understandable before clever

- 📝🔎 Where a scheme lives and what it does is documented in `code_analysis/`, next to the code.
- 🧾📌 Every interface documents its switches, defaults and reference versions in its own `README.md`, and its open work in a `TODO.md`.

## 5. 🟣📖 The Julia twin is a reader's edition

- 🟣🪞 `julia/` holds a hand-written Julia twin of `micro_pumas_v1` and the parts of `micro_pumas_utils` it uses, written for reading and experimenting offline.
- ⚖️🧪 It is verified numerically against the Fortran (standalone column driver in `tests/column/`), never run inside CESM; the Fortran stays the code CESM runs.
- 🔗🏷️ Every Julia function names its Fortran twin in its docstring, and `julia/TWINS.csv` lists the pairs, so a change on one side can be carried to the other.
- 🧮🔁 Same arithmetic in the same order as the Fortran: no physics changes, no fused multiply-add; oddities found while translating go to `julia/FINDINGS.md`.
