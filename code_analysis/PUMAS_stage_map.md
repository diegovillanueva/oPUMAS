# 🗺️🧊 Stage map of `micro_pumas_tend`

📍📂 A reading guide to the main routine of `micro_pumas_v1.F90` (branch `dev`, PUMAS names, as of 2026-10-01).
🎯📖 It replaces a structural refactor: the code stays identical to upstream, and this map shows where each stage lives.
🔎⚓ Line numbers drift with upstream merges, so every stage also has a search anchor, a comment that exists verbatim in the code.
🧮📌 The routine has no sub-step loop; the host calls it once per sub-step (`micro_mg_num_steps`).

## 🪜📋 The stages in order

| # | Stage | Lines | Search anchor | What happens |
|---|---|---|---|---|
| 0 | Interface and declarations | 578–1300 | `! input arguments` | Arguments (inputs from 681, outputs from 724) and about 280 local work arrays (from 821). |
| 1 | Setup | 1305–1710 | `! Initialize local variables` | Local copies of the state, air-density correction of the fall-speed parameters, saturation (`qsat_water`, `qsat_ice`), zeroing of every output and process rate. |
| 2 | Activation and ice nucleation | 1715–1783 | `! calculate potential for droplet activation` | Droplet activation from the host's `npccn` is added to `nc`. Ice nucleation (`nnuccd`, `mnuccd`) when activated nuclei exist below -5 C. |
| 3 | Instantaneous melting and freezing, first pass | 1784–1898 | `! calculate instantaneous precip processes` | Snow and graupel melt at +2 C, rain freezes at -5 C, each limited so the temperature does not cross the threshold. The state is changed directly here. |
| 4 | In-cloud values and fractions | 1899–2027 | `! for microphysical process calculations` | Constant-number switches (`nccons`, `nicons`), precipitation fraction, conversion from grid-mean to in-cloud and in-precipitation values (`qcic`, `qric`, ...), size distributions of liquid and rain. |
| 5 | Process rates, warm rain | 2028–2165 | `! autoconversion of cloud liquid water to rain` | Autoconversion by `kk2000`, the TAU bin scheme, its emulator, or `sb2001`, chosen by `micro_pumas_warm_rain`. |
| 6 | Process rates, ice, snow and mixed phase | 2166–2530 | `! Get size distribution parameters for cloud ice` | Ice autoconversion, size distributions of ice, snow and graupel or hail, immersion and contact freezing, snow aggregation, riming, secondary ice, rain freezing, accretion, Bergeron process, vapour deposition and sublimation. |
| 7 | Process rates, graupel or hail | 2531–2583 | `! Process rate calls for graupel` | Collection of snow, cloud water and rain by graupel, rime splintering, and evaporation and sublimation of precipitation (with or without graupel). |
| 8 | Conservation limiters | 2584–2969 | `! conservation to ensure no negative values` | Every sink of each species is scaled by one ratio so the species cannot go negative in the time step: `qc`, vapour, `nc`, rain, `ni`, snow, graupel. The saturation limit on evaporation and deposition follows from 2906. |
| 9 | Tendency assembly | 2970–3172 | `! get tendencies due to microphysical conversion processes` | The limited rates are summed into `qvlat`, `tlat` and the mass and number tendencies, multiplied by cloud or precipitation fraction. Process rates are stored in `proc_rates` for history. |
| 10 | Re-apply early adjustments | 3173–3274 | `! Assign variables back to start-of-timestep values` | Stages 2 and 3 changed the state directly; here those changes are moved into the tendencies, and the state is updated for sedimentation. |
| 11 | Fall speeds | 3275–3510 | `! obtain new slope parameter to avoid possible singularity` | Mass- and number-weighted fall speeds of all species, with the IFS option and the fall-speed correction (`precip_fall_corr`). |
| 12 | Sedimentation | 3511–3728 | `! Sedimentation` | Implicit (`Sedimentation_implicit`, from 3517) or explicit sub-stepped (`Sedimentation`, from 3620), chosen by `micro_mg_implicit_fall`. Surface precipitation fluxes are summed from 3709. |
| 13 | Instantaneous melting and freezing, second pass | 3729–4022 | `! get new update for variables that includes sedimentation tendency` | Melting of snow and graupel, freezing of rain, melting of cloud ice, homogeneous freezing of droplets at -40 C, and removal of excess supersaturation, all on the state after sedimentation. |
| 14 | Effective radii and final size distributions | 4023–4279 | `! calculate effective radius for pass to radiation code` | Effective radii of droplets and ice for radiation, final size distributions of rain, snow and graupel. |
| 15 | Output diagnostics | 4280–4547 | `! DO STUFF FOR OUTPUT:` | Mean diameters, CloudSat-like and 10 cm radar reflectivity, ice fraction. |

## 🧭🔑 Five things that make the routine readable

- 🔁🧊 **Stages 5 to 7 only compute rates.** Nothing is applied to the state there. Each call fills one or two rate arrays such as `prc` or `psacws`. The physics of each rate lives in `micro_pumas_utils.F90`, in a subroutine with the same name as the call.
- ⚖️✂️ **Stage 8 is where rates stop being independent.** If the sinks of a species exceed what is available, all its sinks are scaled by the same ratio. A new process must be added to the matching limiter block, otherwise mass is not conserved.
- ➕📊 **Stage 9 is the single place where rates become tendencies.** A new process needs one line in each affected budget there (mass, number, latent heat).
- ↩️🌡️ **Stages 3 and 10 are a pair.** The early melting and freezing changes the working state so the rates see it, and stage 10 moves it into the tendencies. Stage 13 repeats the same physics after sedimentation.
- 🧊🔀 **Graupel and hail share one set of variables** (`qg`, `ng`, `lamg`). `do_hail` only selects the particle properties (`mg_hail_props` or `mg_graupel_props`) and density.

## ➕🧩 Checklist for adding a process

1. 🧮 Write the rate as a subroutine in `micro_pumas_utils.F90`, modelled on a neighbour such as `immersion_freezing`.
2. 📞 Call it in stage 5, 6 or 7, behind its switch.
3. ⚖️ Add it to the conservation limiter of every species it removes (stage 8).
4. ➕ Add it to the tendencies and to the latent heating (stage 9), and to `proc_rates` if it should be written.
5. 🎚️ Wire the switch through `micro_pumas_init` and the host interface, off by default, and add a row to `tests/switches.csv`.

## 🔗📚 Related notes

- ❄️ Immersion freezing in detail: `PUMAS_immersion_freezing.md`.
- 🧪 How to check a change: `../tests/README.md`.
