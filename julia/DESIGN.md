# 🧊🪜 Container design for the time step

🙋📝 Written as the M4 draft for Diego's review; the proposed answers of section 4 were adopted. Section 5 records how it was built.
🧭📌 Stage numbers refer to `code_analysis/PUMAS_stage_map.md`.

## 1. 🧱📐 Principles

- 📏 **Struct of arrays.** Every field is a `Matrix{Float64}` of size `(column, level)`, level 1 at the top, like the Fortran. Interface fields are `Vector` (surface) or `(column, level + 1)` (fluxes). The golden files then compare field by field, with no reshaping.
- 🔬 **Scalar physics, looping stages.** The M3 process functions stay scalar. A stage is a plain function `stage_name!(work, parameters)` that loops `for k in levels, i in columns` (Fortran order) and stores the scalar results.
- 🗃️ **Species first, then quantity.** Paths read like the plan's style contract: `state.rain.q`, `incloud.snow.λ`, `rates.rain.accretion_by_snow.mass`.
- 🧊 **Immutable structs of mutable arrays.** Allocated once by `allocate(columns, levels)`, then overwritten in place every step. No module globals; the about 280 Fortran work arrays become fields or stage locals.
- 🔤 **Names from the dictionary.** Each rate's path follows the CF-CCPP name in `translation/dictionary.csv` (`local_rain_mass_tend_accretion_by_snow` → `rates.rain.accretion_by_snow.mass`). A dictionary `Julia` column records the path (plan, rule 9, done on `main_cam_interface`).

## 2. 📦🧩 The seven containers

| Container | Fortran counterpart | Written by stages | Read by stages |
|---|---|---|---|
| `ColumnInputs` | `intent(in)` arguments of `micro_pumas_tend` | host | 1–15 |
| `State` | working copies `t`, `q`, `qc`, `nc`, … | 1, 3, 10, 12, 13 | all |
| `InCloudValues` | `qcic`, `qric`, `lamr`, `n0r`, `umr`, `unr`, fractions, … | 4, 6, 11 | 5–8, 11–14 |
| `LocalRates` | `prc`, `pra`, `psacws`, … (in-cloud or in-precipitation) | 2, 5–8 | 8, 9 |
| `ProcessRates` | `proc_rates_type` (grid mean, for history) | 9, 12, 13 | host, tests |
| `Tendencies` | `tlat`, `qvlat`, `qctend`, …, `*sedten` | 9, 10, 12, 13 | host |
| `Diagnostics` | `effc`, `refl`, `prect`, fluxes, `*out` | 12, 14, 15 | host |

🧺 A `Workspace` bundles all seven, so a stage takes `(work, parameters)` and unpacks only what it needs: `(; state, incloud) = work`.

## 3. 🔍🧾 Field sketches

```julia
"Mass q (kg kg⁻¹) and number N (kg⁻¹) of one hydrometeor."
struct Species
    q::Matrix{Float64}
    N::Matrix{Float64}
end

"Inputs from the host, never modified (`intent(in)` of micro_pumas_tend)."
struct ColumnInputs
    time_step::Float64                        # deltatin
    T::Matrix{Float64}                        # t
    q_vapor::Matrix{Float64}                  # q
    cloud_liquid::Species                     # qcn, ncn
    cloud_ice::Species                        # qin, nin
    rain::Species; snow::Species; graupel::Species
    pressure::Matrix{Float64}                 # p
    pressure_thickness::Matrix{Float64}       # pdel
    interface_pressure::Matrix{Float64}       # pint, levels + 1
    cloud_fraction::Matrix{Float64}           # cldn
    liquid_cloud_fraction::Matrix{Float64}    # liqcldf
    ice_cloud_fraction::Matrix{Float64}       # icecldf
    saturation_scaling::Matrix{Float64}       # qsatfac
    relative_variance::Matrix{Float64}        # relvar
    accretion_enhancement::Matrix{Float64}    # accre_enhan
    activated_droplets::Matrix{Float64}       # npccn
    activated_ice_nuclei::Matrix{Float64}     # naai
    host_freezing::NamedTuple                 # frzimm, frzcnt, frzdep
    snow_from_host::Species                   # tnd_qsnow, tnd_nsnow (as rates)
    ice_effective_radius::Matrix{Float64}     # re_ice
end

"Working copy of the state, changed by stages 3, 10, 12 and 13."
struct State
    T::Matrix{Float64}
    q_vapor::Matrix{Float64}
    cloud_liquid::Species; cloud_ice::Species
    rain::Species; snow::Species; graupel::Species
    air_density::Matrix{Float64}              # rho
    q_sat_liquid::Matrix{Float64}             # qvl
    q_sat_ice::Matrix{Float64}                # qvi
    air::NamedTuple   # diffusivity, viscosity, schmidt, fall speed prefactors a
end

"In-cloud or in-precipitation values of one species plus its size distribution."
struct InCloudSpecies
    q::Matrix{Float64}; N::Matrix{Float64}    # qric, nric
    λ::Matrix{Float64}; N₀::Matrix{Float64}   # lamr, n0r
    fall_speed_mass::Matrix{Float64}          # umr
    fall_speed_number::Matrix{Float64}        # unr
end
# 💧 cloud liquid adds μ (pgam); cloud ice and liquid have no N₀ at this stage.

struct InCloudValues
    cloud_liquid::InCloudSpecies; cloud_ice::InCloudSpecies
    rain::InCloudSpecies; snow::InCloudSpecies; graupel::InCloudSpecies
    precipitation_fraction::Matrix{Float64}   # precip_frac
    liquid_fraction::Matrix{Float64}          # lcldm
    ice_fraction::Matrix{Float64}             # icldm
end

"One local rate pair (mass kg kg⁻¹ s⁻¹, number kg⁻¹ s⁻¹)."
struct Rate
    mass::Matrix{Float64}
    number::Matrix{Float64}
end

# 🧮 LocalRates: species that loses the mass → process → Rate, following the
#    CF-CCPP dictionary names. Built as nested NamedTuples by allocate(), e.g.
#      rates.cloud_liquid.autoconversion       (prc, nprc)
#      rates.cloud_liquid.accretion_by_rain     (pra, npra)
#      rates.cloud_liquid.accretion_by_snow     (psacws, npsacws)
#      rates.rain.accretion_by_snow             (pracs, npracs)
#      rates.rain.freezing                      (mnuccr, nnuccr)
#      rates.cloud_ice.autoconversion           (prci, nprci)
#      rates.snow.self_aggregation.number       (nsagg)
#      rates.vapor.deposition_onto_snow.mass    (vap_deps)
#    with `number` absent (or a shared zero) for mass-only processes.

"proc_rates_type: grid-mean process rates kept for history, same field order."
struct ProcessRates   # 93 fields, from the dictionary rows of proc_rates_type
    ...
end

struct Tendencies
    heating::Matrix{Float64}                  # tlat (W kg⁻¹)
    q_vapor::Matrix{Float64}                  # qvlat
    cloud_liquid::Species; cloud_ice::Species # qctend, nctend, ...
    rain::Species; snow::Species; graupel::Species
    sedimentation::NamedTuple                 # qcsedten ... ngsedten as Species
end

struct Diagnostics   # effc, effi, reff_*, prect, preci, *flx, refl..., *out, freq*
    ...
end
```

## 4. 🙋❓ Questions for Diego

1. 🗃️ **Rates by losing species or by process?** Proposed: `rates.rain.accretion_by_snow.mass`, which matches the dictionary. The alternative is flat process names such as `rates.rain_accretion_by_snow.mass`, which is shorter to unpack but drifts from the dictionary.
2. 📦 **Seven containers in one `Workspace`?** Proposed yes, so every stage has the same signature. The alternative passes only the needed containers, which documents data flow in each signature but makes `step.jl` noisier.
3. 🧾 **`ProcessRates` keeps PUMAS order and meaning?** Proposed: same fields as `proc_rates_type` with dictionary names, so history output and golden comparison stay one-to-one. Its 93 fields would be generated from `translation/dictionary.csv` once and then written out by hand-readable code, not by macros.
4. 🔣 **Unicode field names?** Proposed: `λ`, `N₀`, `μ` only inside `InCloudSpecies`, where they are standard; everything else is words. The alternative is `slope`, `intercept`, `shape` everywhere, as in `size_distribution`'s return value today.
5. 🧭 **`air` as a NamedTuple of matrices in `State`?** It holds `diffusivity`, `viscosity`, `schmidt` and the density-corrected fall speed prefactors. The vapor twins of M3 already take a scalar `air` NamedTuple, built per point in the stage loop.

## 5. ✅🧱 As built (M4–M6)

- 📦 **NamedTuples, not structs.** The containers are nested NamedTuples of matrices built by `allocate` from the path tables in `src/names.jl` (`LOCAL_RATE_PATHS`, `PROCESS_RATE_PATHS`, `TENDENCY_PATHS`, `DIAGNOSTIC_PATHS`, `FLUX_PATHS`, `SURFACE_PATHS`), so a new field is one table row. `work` is `(; inputs, state, incloud, rates, process_rates, tendencies, diagnostics, scratch)`; `scratch` holds the Fortran work arrays that several stages share (`dumc`, `fc`, `dum_2D`, …).
- 🔤 **Readable locals in the stages (M6).** Each stage starts with `(; inputs, state, ...) = work` and binds the fields it uses to readable locals with the PUMAS name in a comment, `qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc`. Long paths go through a species sub-container (`liquid_means = process_rates.cloud_liquid`). The M4 `work.fortran` aliases with PUMAS names are gone; `SYMBOLS.md` lists every binding.
- 🔣 **Constants under their equation symbols.** `host_constants(parameters)` gives `g`, `Rd`, `cp`, `Lv`, `Lf`, `Ls`, `T₀`, …; the PUMAS short constants (`qsmall`, `omsm`, `mincld`) are the long names of `src/constants.jl`.
- ⚖️ **Numbers first, names second.** Every renaming commit kept all 16 switch rows bit-for-bit; renaming and arithmetic never shared a commit.
