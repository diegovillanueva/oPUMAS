"""
🟣📖 oPUMAS: a readable Julia twin of the PUMAS cloud microphysics (micro_pumas_v1.F90).

🧊🔗 Same arithmetic, same order as the Fortran; names follow the CF/CCPP vocabulary (translation/dictionary.csv).
"""
module oPUMAS

using NCDatasets

include("libm.jl")
include("gamma.jl")
include("constants.jl")
include("saturation.jl")
include("hydrometeors.jl")
include("parameters.jl")
include("namelist.jl")
include("size_distributions.jl")
include("processes/warm_rain.jl")
include("processes/bin_warm_rain.jl")
include("processes/aggregation.jl")
include("processes/riming.jl")
include("processes/freezing.jl")
include("processes/deposition.jl")
include("processes/graupel.jl")
include("processes/effective_radius.jl")
include("names.jl")
include("workspace.jl")
include("stages/common.jl")
include("stages/s01_setup.jl")
include("stages/s02_activation.jl")
include("stages/s03_instant_phase_change.jl")
include("stages/s04_incloud.jl")
include("stages/s05_warm_rain.jl")
include("stages/s06_ice_processes.jl")
include("stages/s07_graupel_processes.jl")
include("stages/s08_conservation.jl")
include("stages/s09_tendencies.jl")
include("stages/s10_reapply.jl")
include("stages/s11_fall_speeds.jl")
include("stages/s12_sedimentation.jl")
include("stages/s13_instant_phase_change2.jl")
include("stages/s14_effective_radius.jl")
include("stages/s15_output_diagnostics.jl")
include("step.jl")
include("time_loop.jl")

export gamma_function
export SaturationConstants, saturation_over_water, saturation_over_ice
export HydrometeorProperties, hydrometeor_properties
export HostConstants, Settings, Parameters, setup, read_settings
export size_distribution, droplet_size_distribution, mean_diameter
export rising_factorial, variance_coefficient, latent_heat_factor
export kk2000_autoconversion, sb2001_autoconversion, sb2001_accretion, accretion
export rain_self_collection
export tau_bins, tau_bin_distribution!, tau_collision_rates!, tau_collection
export ice_to_snow_autoconversion, snow_self_aggregation, ice_collection_by_snow
export droplet_collection_by_snow, rime_splintering_by_snow, rain_collection_by_snow
export rain_immersion_freezing
export ventilated_diffusion, bergeron_to_snow, snow_deposition, ice_deposition
export precipitation_evaporation
export snow_collection_by_rain, droplet_collection_by_graupel, rimed_snow_to_graupel
export rain_collection_by_graupel, rain_snow_collisions_to_graupel
export rime_splintering_by_graupel, add_rain_cloud_effective_radius
export ColumnInputs, read_inputs, allocate, output_paths, field_at
export microphysics_step!, advance_state!, time_loop!

end
