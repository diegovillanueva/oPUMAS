"""
🟣📖 oPUMAS: a readable Julia twin of the PUMAS cloud microphysics (micro_pumas_v1.F90).

🧊🔗 Same arithmetic, same order as the Fortran; names follow the CF-CCPP vocabulary.
"""
module oPUMAS

include("gamma.jl")
include("constants.jl")
include("saturation.jl")
include("hydrometeors.jl")
include("parameters.jl")
include("namelist.jl")
include("size_distributions.jl")
include("processes/warm_rain.jl")

export gamma_function
export SaturationConstants, saturation_over_water, saturation_over_ice
export HydrometeorProperties, hydrometeor_properties
export HostConstants, Settings, Parameters, setup, read_settings
export size_distribution, droplet_size_distribution, mean_diameter
export rising_factorial, variance_coefficient, latent_heat_factor
export kk2000_autoconversion, sb2001_autoconversion, sb2001_accretion, accretion
export rain_self_collection

end
