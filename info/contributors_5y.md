# 👥🧊 ESCOMP/PUMAS contributors, last 5 years

🗓️📊 Window: PRs created 2021-09-30 → 2026-09-30 (`gh pr list --repo ESCOMP/PUMAS --search "created:>=2021-09-30"`), 32 PRs opened, 24 merged.
🔎🧾 Identities come from public GitHub profiles and upstream commit author emails (`git log --all --since=2021-09-30 --format='%an <%ae>'`); anything marked *inferred* is not stated on the profile.

## 📈🏆 Ranking

| # | GitHub | Person | Affiliation | PRs opened | Merged | PRs merged by them |
|---|--------|--------|-------------|-----------:|-------:|-------------------:|
| 1 | Katetc | Kate Thayer-Calder | NCAR, Boulder | 14 | 10 | 21 |
| 2 | cacraigucar | Cheryl Craig (from commit author name) | NCAR/UCAR (`cacraig@ucar.edu`) | 7 | 6 | 2 |
| 3 | sjsprecious | Jian Sun | NCAR, Boulder | 6 | 6 | 1 |
| 4 | nusbaume | Jesse Nusbaumer | NCAR CGD, AMP section | 1 | 1 | 0 |
| 4 | andrewgettelman | Andrew Gettelman | PNNL (formerly NCAR) | 1 | 1 | 0 |
| – | djgagne | David John Gagne | NSF NCAR, head of MILES group | 1 | 0 (open) | 0 |
| – | jimmielin | Haipeng Lin | NCAR (`hplin@ucar.edu`) | 1 | 0 (open) | 0 |
| – | wkchuang | W. Chuang (*inferred*: Columbia, `wc2227@columbia.edu`, NCAR HPC user) | Columbia University | 1 | 0 (closed) | 0 |

## 🧑‍🔬📋 Who they are and what they did

### 🔑👑 Kate Thayer-Calder (Katetc): maintainer and gatekeeper
- 🏢🌨️ NCAR scientist working on atmosphere and land-ice models in CESM.
- 🧩🔀 Merged 21 of the 24 PRs, so she is the reviewer any oPUMAS PR must convince.
- 🛠️📦 Her PRs: PUMAS v1 renaming (#37), het-freezing updates (#48), ML warm-rain tau emulator port and output fields (#53, #57–#60), vapor deposition onto snow (#63), README (#69), bugfix v1.39 (#74).

### 🧰🔗 Cheryl Craig (cacraigucar): CAM software engineering and integration
- 🏢💻 NCAR/UCAR CAM software engineer (profile empty; name from commit author, and the author of the most upstream commits since 2021 (77)).
- 🔄🧱 Her PRs: number-concentration update (#41), nt2/Bigg v2 (#42), PUMAS DDT (#46), extra ML mods (#54), CCPP round 3 (#75), `_r8` fixes (#78, the current upstream HEAD).

### ⚡🖥️ Jian Sun (sjsprecious): GPU/OpenACC
- 🏢🧮 NCAR software engineer: HPC, GPU programming, profiling, optimisation.
- 🚀🐞 All 6 PRs keep PUMAS working on GPUs: OpenACC port and fixes (#34, #35, #39, #43, #56, #68).

### 🧱🔌 Jesse Nusbaumer (nusbaume): CCPP-ization
- 🏢🌍 NCAR CGD software engineer (AMP section).
- 🧩📐 #71 removed CAM dependencies and added the CCPP layer for PUMAS v1, the reason for `micro_pumas_ccpp.F90/.meta`.

### 🌧️📚 Andrew Gettelman (andrewgettelman): scientific originator
- 🏢🔬 Senior Earth Scientist at PNNL, previously NCAR; the "G" in Morrison–Gettelman (MG) microphysics that PUMAS descends from.
- 📡🧊 #36 changed the radar reflectivity calculation.

### 🤖🧠 Machine-learning contributors (open or closed PRs)
- 🧭🔬 David John Gagne (djgagne): heads NSF NCAR's MILES ML group; #73 "1nn tau emulator updates" is still open.
- 🧪📉 W. Chuang (wkchuang, *inferred* Columbia): #65 made the tau emulator output tendencies (closed unmerged).

### 🏷️📑 Haipeng Lin (jimmielin): CCPP metadata
- 🏢🧾 NCAR (commit email `hplin@ucar.edu`); #76 updates the het-freezing standard name in `micro_pumas_ccpp.meta` (open).

## 🧭💡 Takeaways for an oPUMAS PR
- 🎯👩‍🔬 Plan around Kate Thayer-Calder as reviewer; she merges almost everything.
- 🧱🔌 Code must respect the CCPP interface (`micro_pumas_ccpp.*`, standard names), see #71/#75/#76.
- ⚡🖥️ Changes inside kernels need to stay OpenACC-safe; Jian Sun's PRs show GPU breakage is a recurring theme.
- 🐢📉 About 6 PRs a year, all from NCAR staff so far: an external PR is unusual, so contact the maintainers (issue or email) before submitting.
