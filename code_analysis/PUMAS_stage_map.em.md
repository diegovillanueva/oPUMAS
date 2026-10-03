# 🗺️🧊 Stage map of `micro_pumas_tend` — 🖥️👓 Vim edition

📜🔁 Plain-text twin of `PUMAS_stage_map.md`: no tables, no wide lines, readable in Vim.
📍📂 Routine: `micro_pumas_tend` in `micro_pumas_v1.F90` (branch `main`, then called `main_cam_interface`, 2026-10-01).
🎯📖 The code stays identical to upstream; this map shows where each stage lives.
🔎⚓ Line numbers drift with merges, so jump with the anchor: `/! Sedimentation` in Vim.
🧮📌 No sub-step loop inside: the host calls it once per sub-step (`micro_mg_num_steps`).

🔤🎨 Legend
    🟦 setup / bookkeeping      🟩 rates only (state untouched)
    🟥 state changed directly   🟨 rates → tendencies         🟪 diagnostics


## 🌊🧭 Flow at a glance

```
  🟦  0  📜 Interface & declarations ........................  578
  🟦  1  🧰 Setup: copies, qsat, zero everything ............. 1305
  🟥  2  💧 Activation & ice nucleation ...................... 1715
  🟥  3  🔥❄️ Instant melt/freeze, 1st pass ──────────────┐ .. 1784
  🟦  4  ☁️ In-cloud values & fractions                   │ .. 1899
  🟩  5  🌧️ Warm rain rates                               │ .. 2028
  🟩  6  ❄️ Ice / snow / mixed-phase rates                │ .. 2166
  🟩  7  🧊 Graupel / hail rates                          │ .. 2531
  ⚖️  8  ✂️ Conservation limiters                         │ .. 2584
  🟨  9  ➕ Tendency assembly                             │ .. 2970
  🟨 10  ↩️ Re-apply early adjustments  ◄─────────────────┘ .. 3173
  🟦 11  🪂 Fall speeds ...................................... 3275
  🟥 12  ⬇️ Sedimentation .................................... 3511
  🟥 13  🔥❄️ Instant melt/freeze, 2nd pass (post-sedim.) .... 3729
  🟪 14  📏 Effective radii & final size distributions ....... 4023
  🟪 15  📡 Output diagnostics ............................... 4280
                                              🏁 end subroutine 4547
```


## 🪜📋 The stages, one by one

### 0️⃣ 📜🟦 Interface and declarations
    📍 lines   578–1300
    ⚓ anchor  `! input arguments`
    📝 what    🔌 Arguments: inputs from 681, outputs from 724.
               🗃️ About 280 local work arrays from 821.

### 1️⃣ 🧰🟦 Setup
    📍 lines   1305–1710
    ⚓ anchor  `! Initialize local variables`
    📝 what    📋 Local copies of the state.
               🌬️ Air-density correction of the fall-speed parameters.
               💦 Saturation: `qsat_water`, `qsat_ice`.
               0️⃣ Zeroing of every output and process rate.

### 2️⃣ 💧🟥 Activation and ice nucleation
    📍 lines   1715–1783
    ⚓ anchor  `! calculate potential for droplet activation`
    📝 what    💧 Droplet activation from the host's `npccn` is added to `nc`.
               ❄️ Ice nucleation (`nnuccd`, `mnuccd`) when activated nuclei exist below -5 C.

### 3️⃣ 🔥❄️🟥 Instantaneous melting and freezing, first pass
    📍 lines   1784–1898
    ⚓ anchor  `! calculate instantaneous precip processes`
    📝 what    🫠 Snow and graupel melt at +2 C.
               🧊 Rain freezes at -5 C.
               🚧 Each is limited so the temperature does not cross the threshold.
               ⚠️ The state is changed directly here (undone into tendencies in stage 10).

### 4️⃣ ☁️🟦 In-cloud values and fractions
    📍 lines   1899–2027
    ⚓ anchor  `! for microphysical process calculations`
    📝 what    🔒 Constant-number switches (`nccons`, `nicons`).
               🌂 Precipitation fraction.
               🔄 Grid-mean → in-cloud / in-precip values (`qcic`, `qric`, ...).
               📐 Size distributions of liquid and rain.

### 5️⃣ 🌧️🟩 Process rates, warm rain
    📍 lines   2028–2165
    ⚓ anchor  `! autoconversion of cloud liquid water to rain`
    📝 what    🔀 Autoconversion chosen by `micro_pumas_warm_rain`:
                  • `kk2000`
                  • TAU bin scheme
                  • TAU emulator
                  • `sb2001`

### 6️⃣ ❄️🟩 Process rates, ice, snow and mixed phase
    📍 lines   2166–2530
    ⚓ anchor  `! Get size distribution parameters for cloud ice`
    📝 what    🧊 Ice autoconversion.
               📐 Size distributions of ice, snow and graupel or hail.
               🥶 Immersion and contact freezing.
               ❄️ Snow aggregation, riming, secondary ice.
               🌧️ Rain freezing, accretion.
               🔁 Bergeron process.
               💨 Vapour deposition and sublimation.

### 7️⃣ 🧊🟩 Process rates, graupel or hail
    📍 lines   2531–2583
    ⚓ anchor  `! Process rate calls for graupel`
    📝 what    🧲 Collection of snow, cloud water and rain by graupel.
               ✨ Rime splintering.
               💨 Evaporation and sublimation of precipitation (with or without graupel).

### 8️⃣ ⚖️✂️ Conservation limiters
    📍 lines   2584–2969
    ⚓ anchor  `! conservation to ensure no negative values`
    📝 what    ✂️ Every sink of a species is scaled by one ratio so it cannot go negative.
               🔢 Order: `qc` → vapour → `nc` → rain → `ni` → snow → graupel.
               💦 Saturation limit on evaporation and deposition from 2906.

### 9️⃣ ➕🟨 Tendency assembly
    📍 lines   2970–3172
    ⚓ anchor  `! get tendencies due to microphysical conversion processes`
    📝 what    ➕ Limited rates summed into `qvlat`, `tlat`, mass and number tendencies.
               ☁️ Each multiplied by cloud or precipitation fraction.
               🗄️ Rates stored in `proc_rates` for history output.

### 🔟 ↩️🟨 Re-apply early adjustments
    📍 lines   3173–3274
    ⚓ anchor  `! Assign variables back to start-of-timestep values`
    📝 what    ↩️ Stages 2 and 3 changed the state; those changes move into the tendencies.
               🔄 State is updated for sedimentation.

### 1️⃣1️⃣ 🪂🟦 Fall speeds
    📍 lines   3275–3510
    ⚓ anchor  `! obtain new slope parameter to avoid possible singularity`
    📝 what    🪂 Mass- and number-weighted fall speeds of all species.
               🎛️ IFS option and fall-speed correction (`precip_fall_corr`).

### 1️⃣2️⃣ ⬇️🟥 Sedimentation
    📍 lines   3511–3728
    ⚓ anchor  `! Sedimentation`
    📝 what    🔀 Chosen by `micro_mg_implicit_fall`:
                  • implicit `Sedimentation_implicit` (from 3517)
                  • explicit sub-stepped `Sedimentation` (from 3620)
               🌍 Surface precipitation fluxes summed from 3709.

### 1️⃣3️⃣ 🔥❄️🟥 Instantaneous melting and freezing, second pass
    📍 lines   3729–4022
    ⚓ anchor  `! get new update for variables that includes sedimentation tendency`
    📝 what    🫠 Melting of snow and graupel, melting of cloud ice.
               🧊 Freezing of rain; homogeneous freezing of droplets at -40 C.
               💨 Removal of excess supersaturation.
               🔁 All on the state after sedimentation.

### 1️⃣4️⃣ 📏🟪 Effective radii and final size distributions
    📍 lines   4023–4279
    ⚓ anchor  `! calculate effective radius for pass to radiation code`
    📝 what    ☀️ Effective radii of droplets and ice for radiation.
               📐 Final size distributions of rain, snow and graupel.

### 1️⃣5️⃣ 📡🟪 Output diagnostics
    📍 lines   4280–4547
    ⚓ anchor  `! DO STUFF FOR OUTPUT:`
    📝 what    📏 Mean diameters.
               📡 CloudSat-like and 10 cm radar reflectivity.
               🧊 Ice fraction.


## 🧭🔑 Five things that make the routine readable

🔁🟩 **Stages 5–7 only compute rates.**
    🚫 Nothing is applied to the state there.
    🗂️ Each call fills one or two rate arrays such as `prc` or `psacws`.
    📚 The physics lives in `micro_pumas_utils.F90`, in a subroutine named like the call.

⚖️✂️ **Stage 8 is where rates stop being independent.**
    📉 If the sinks of a species exceed what is available, all are scaled by one ratio.
    ⚠️ A new process must join the matching limiter block, or mass is not conserved.

➕📊 **Stage 9 is the single place where rates become tendencies.**
    ✍️ A new process needs one line per affected budget: mass, number, latent heat.

↩️🌡️ **Stages 3 and 10 are a pair.**
    🔥 Early melt/freeze changes the working state so the rates see it.
    ↩️ Stage 10 moves it into the tendencies.
    🔁 Stage 13 repeats the same physics after sedimentation.

🧊🔀 **Graupel and hail share one set of variables** (`qg`, `ng`, `lamg`).
    🎚️ `do_hail` only selects particle properties (`mg_hail_props` / `mg_graupel_props`)
       and density.


## ➕🧩 Checklist for adding a process

    [ ] 1️⃣ 🧮 Write the rate as a subroutine in `micro_pumas_utils.F90`
              (model it on a neighbour such as `immersion_freezing`).
    [ ] 2️⃣ 📞 Call it in stage 5, 6 or 7, behind its switch.
    [ ] 3️⃣ ⚖️ Add it to the limiter of every species it removes (stage 8).
    [ ] 4️⃣ ➕ Add it to the tendencies and latent heating (stage 9),
              and to `proc_rates` if it should be written.
    [ ] 5️⃣ 🎚️ Wire the switch through `micro_pumas_init` and the host interface,
              off by default, and add a row to `tests/switches.csv`.


## 🔗📚 Related notes

    ❄️🔬 Immersion freezing in detail   →  `PUMAS_immersion_freezing.md`
    🧪✅ How to check a change          →  `../tests/README.md`
    📄📊 Table version of this map      →  `PUMAS_stage_map.md`
