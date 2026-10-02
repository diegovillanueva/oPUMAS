# 🏷️🔣 PUMAS short names of the constants used by the M4 stages, so that each stage reads
#    like its Fortran block. The M6 readability pass replaces them by the long names.

const qsmall = smallest_mixing_ratio
const omsm = one_minus_small
const mincld = minimum_cloud_fraction
const mi0 = new_ice_crystal_mass
const mi0l_min = minimum_frozen_droplet_mass
const rhow = bulk_density.liquid_water
const rhoi = bulk_density.cloud_ice
const rhows = bulk_density.solid_water
const rhosn = bulk_density.snow
const sublim_factor = number_sublimation_factor

"Columns and levels of a workspace (`mgncol`, `nlev`)."
grid_size(work) = size(work.inputs.T)

"""
    host_short_names(parameters)

🌍🔣 Host constants and derived temperatures under their PUMAS names.
"""
function host_short_names(parameters)
    h = parameters.host
    return (g = h.gravity, r = h.gas_constant_dry_air, rv = h.gas_constant_water_vapor,
        cpp = h.specific_heat_dry_air, tmelt = h.melting_point,
        xxlv = h.latent_heat_vaporization, xlf = h.latent_heat_fusion,
        xxls = parameters.latent_heat_sublimation, rhosu = parameters.reference_air_density,
        icenuct = parameters.ice_nucleation_temperature,
        snowmelt = parameters.snow_melting_temperature,
        rainfrze = parameters.rain_freezing_temperature)
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
