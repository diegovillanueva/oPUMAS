# 🧭🧊 oPUMAS philosophy: what separates oPUMAS from PUMAS

🌱📘 oPUMAS is a development copy of [ESCOMP/PUMAS](https://github.com/ESCOMP/PUMAS) that aims to be easier to read, learn and extend, without drifting away from PUMAS.
🧩📐 Every design decision in this repo is checked against the principles below; a new principle is added here before it is applied in the code.

## 1. 🔤📖 Names are for readers, the dictionary is for the sync

- 🔤✨ oPUMAS uses variable and function names that minimise the effort of understanding what they mean.
- 📖🔁 Every renamed identifier keeps a dictionary entry that maps it to its reference name in PUMAS, so that syncing PUMAS commits into oPUMAS (and sending changes back) stays easy.
- 📄🗂️ The dictionary is [`name_map.csv`](name_map.csv), one row per renamed identifier, with columns `pumas_name, opumas_name, kind, scope, meaning`.
  - 🏷️ `kind`: `variable`, `function`, `subroutine`, `module`, `namelist` or `type`.
  - 📍 `scope`: file or routine in which the name is valid; `global` for names that are unique across the code.
- 🤖🔄 The dictionary is machine-readable on purpose: an upstream diff can be translated to oPUMAS names (and an oPUMAS diff back to PUMAS names) by a script instead of by hand.
- ⚖️🧪 A rename never changes results: a commit that renames does nothing else, and is verified bit-for-bit against the commit before it.
- 🎛️🤝 Names that users and host models see (namelist variables, history field names, the `micro_pumas_init`/`micro_pumas_tend` argument lists) keep their PUMAS spelling, so existing namelists, scripts and host interfaces keep working.
- 🚧📅 Status: no identifier has been renamed yet; `name_map.csv` holds only its header.

## 2. 🎛️🔕 Everything new is switchable and defaults to PUMAS

- 🎚️✅ Every feature added in oPUMAS or in an `interface/*` is switchable and off by default: an empty namelist makes oPUMAS behave exactly like upstream PUMAS, bit-for-bit where possible.
- 🧭📌 When out-of-the-box CESM–PUMAS versions differ, the CESM 2.3 behaviour is the default; settings of other versions are offered as namelist files (e.g. `interface/cesm215/user_nl_cam.pumas_cesm3like`).

## 3. 🔄🌿 Stay close to upstream

- 🔁📥 Upstream PUMAS commits are merged continuously (`main_cam` mirrors `upstream/main_cam`; see [`git_work.md`](git_work.md)).
- 🧱🔌 Host-model glue lives in `interface/<host>/`, never in the oPUMAS core files, so the core stays comparable with PUMAS.

## 4. 📚🪜 Understandable before clever

- 📝🔎 Where a scheme lives and what it does is documented in `code_analysis/`, next to the code.
- 🧾📌 Every interface documents its switches, defaults and reference versions in its own `README.md`, and its open work in a `TODO.md`.
