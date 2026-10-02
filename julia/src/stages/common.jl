# 💧🧊 Bulk densities of the stages under short readable names (`rhow`, `rhoi`, `rhows`,
#    `rhosn` in the Fortran).

const ρ_water = bulk_density.liquid_water
const ρ_cloud_ice = bulk_density.cloud_ice
const ρ_solid_water = bulk_density.solid_water
const ρ_snow = bulk_density.snow

"Columns and levels of a workspace (`mgncol`, `nlev`)."
grid_size(work) = size(work.inputs.T)

"""
    host_constants(parameters)

🌍🔣 Host constants and derived temperatures under the symbols of the equations: `g`,
`Rd` (`r`), `Rv` (`rv`), `cp` (`cpp`), `T₀` (`tmelt`), `Lv` (`xxlv`), `Lf` (`xlf`), `Ls`
(`xxls`), `ρ_ref` (`rhosu`), `T_ice_nucleation` (`icenuct`), `T_snow_melt` (`snowmelt`)
and `T_rain_freeze` (`rainfrze`).
"""
function host_constants(parameters)
    h = parameters.host
    return (g = h.gravity, Rd = h.gas_constant_dry_air, Rv = h.gas_constant_water_vapor,
        cp = h.specific_heat_dry_air, T₀ = h.melting_point,
        Lv = h.latent_heat_vaporization, Lf = h.latent_heat_fusion,
        Ls = parameters.latent_heat_sublimation, ρ_ref = parameters.reference_air_density,
        T_ice_nucleation = parameters.ice_nucleation_temperature,
        T_snow_melt = parameters.snow_melting_temperature,
        T_rain_freeze = parameters.rain_freezing_temperature)
end

"""
    dense_precipitation(parameters)

🟤🔣 Fall speed law and bulk density of the dense precipitation species (graupel or
hail): `(; a, b, density)` as `ag`/`ah`, `bgtmp` and `rhogtmp` of the Fortran. With
neither species on, `b` and `density` are 0 as in the Fortran.
"""
function dense_precipitation(parameters)
    settings = parameters.settings
    law = settings.do_hail ? fall_speed_law.hail : fall_speed_law.graupel
    density = settings.do_hail ? bulk_density.hail :
              settings.do_graupel ? bulk_density.graupel : 0.0
    return (a = law.a, b = dense_precipitation_exponent(settings), density)
end
