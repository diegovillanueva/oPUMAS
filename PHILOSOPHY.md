# 🧭🧊 oPUMAS philosophy: what separates oPUMAS from PUMAS

🌱📘 oPUMAS is a development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) that aims to be easier to read, learn and extend, without drifting away from PUMAS.
🧩📐 Every design decision in this repo is checked against the principles below; a new principle is added here before it is applied in the code.

## 1. 🔤📖 Names are for readers, the dictionary is for the sync

- 🔤✨ oPUMAS offers variable and function names that minimise the effort of understanding what they mean.
- 📖🔁 Every such name keeps a dictionary entry that maps it to its reference name in PUMAS, so that syncing PUMAS commits into oPUMAS (and sending changes back) stays easy.
- 🧊🔤 The Fortran keeps the upstream PUMAS names. The readable names live in the Julia twin (principle 5), which is the reader's edition of the same code.
- 📄🗂️ The dictionary is [`translation/dictionary.csv`](translation/dictionary.csv): one row per concept, with the `Fortran` name, the `Julia` path, units and meaning. Details: [`translation/README.md`](translation/README.md).
- 🔁✅ Dictionary, Julia twin and Fortran are kept in sync by a test (`julia/test/names_check.jl`, part of `Pkg.test()`): every Julia path exists, every Julia name table entry has its row, and every Fortran name still exists in the sources.
- 🎛️🤝 Names that users and host models see (module names, namelist variables, history field names, `micro_pumas_init`/`micro_pumas_tend` and their argument lists) are never renamed.
- 🗄️🌿 History: until 2026-10-03 a translator kept a second Fortran branch `CF-CCPP` in long names (tag `archive/CF-CCPP`). It was retired because the Julia twin already serves readers, and one Fortran spelling keeps upstream merges simple.
- 🚧📅 Status (2026-10-03): 352 rows (process rates, state, size distributions, fall speeds, sedimentation, outputs of `micro_pumas_tend`, namelist settings).

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
