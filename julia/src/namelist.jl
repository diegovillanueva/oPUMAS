# 📄🔤 Read the `pumas_settings` namelist written by the CAM dump or the column driver.

"CAM namelist name => `Settings` field."
const NAMELIST_FIELDS = Dict(
    "rhmini" => :minimum_rh_for_ice_cloud,
    "rhmini_const" => :minimum_rh_for_ice_cloud,
    "micro_mg_dcs" => :ice_to_snow_size,
    "micro_mg_do_hail" => :do_hail,
    "micro_mg_do_graupel" => :do_graupel,
    "microp_uniform" => :uniform_subcolumns,
    "do_cldice" => :prognostic_cloud_ice,
    "use_hetfrz_classnuc" => :host_heterogeneous_freezing,
    "micro_mg_precip_frac_method" => :precipitation_fraction_method,
    "micro_mg_berg_eff_factor" => :bergeron_efficiency,
    "micro_mg_accre_enhan_fact" => :accretion_enhancement,
    "micro_mg_autocon_fact" => :autoconversion_prefactor,
    "micro_mg_autocon_nd_exp" => :autoconversion_number_exponent,
    "micro_mg_autocon_lwp_exp" => :autoconversion_water_exponent,
    "micro_mg_homog_size" => :homogeneous_freezing_size,
    "micro_mg_vtrmi_factor" => :ice_fall_speed_factor,
    "micro_mg_vtrms_factor" => :snow_fall_speed_factor,
    "micro_mg_effi_factor" => :ice_effective_radius_factor,
    "micro_mg_iaccr_factor" => :ice_accretion_factor,
    "micro_mg_max_nicons" => :maximum_ice_number,
    "allow_sed_supersat" => :remove_supersaturation,
    "micro_pumas_warm_rain" => :warm_rain,
    "micro_pumas_kernel_file" => :tau_kernel_file,
    "micro_mg_evap_sed_off" => :no_sedimentation_evaporation,
    "micro_mg_icenuc_rh_off" => :ice_nucleation_without_rh,
    "micro_mg_icenuc_use_meyers" => :meyers_ice_nucleation,
    "micro_mg_evap_scl_ifs" => :ifs_evaporation_scaling,
    "micro_mg_evap_rhthrsh_ifs" => :ifs_evaporation_threshold,
    "micro_mg_rainfreeze_ifs" => :ifs_rain_freezing,
    "micro_mg_ifs_sed" => :ifs_snow_fall_speed,
    "micro_mg_precip_fall_corr" => :precipitation_fall_correction,
    "micro_mg_accre_sees_auto" => :accretion_sees_autoconversion,
    "micro_mg_implicit_fall" => :implicit_sedimentation,
    "micro_mg_nccons" => :constant_droplet_number,
    "micro_mg_nicons" => :constant_ice_number,
    "micro_mg_ncnst" => :droplet_number,
    "micro_mg_ninst" => :ice_number,
    "micro_mg_ngcons" => :constant_graupel_number,
    "micro_mg_ngnst" => :graupel_number,
    # 🧊🎛️ Sppe controls (micro_pumas_sppe.F90)
    "naai_het_also_in_mpc" => :nucleation_in_mixed_phase,
    "enable_nimax" => :cap_ice_number,
    "detrainment_ramp_liq" => :supersaturation_to_liquid,
    "onoff_contact_freezing" => :droplet_contact_freezing,
    "onoff_sip" => :rime_splintering,
    "rain_freeze" => :rain_freezing,
    "auto_liq_factor" => :autoconversion_factor,
    "riming_factor" => :riming_factor,
    "micro_mg_ice_berg_ice_factor" => :bergeron_ice_factor,
    "micro_mg_ice_berg_snow_factor" => :bergeron_snow_factor,
    "factor_hm" => :splinter_number_factor
)

"Namelist entries that do not change the translated physics."
const NAMELIST_IGNORED = ()

"""
    read_settings(path) -> Settings

Read namelist group `pumas_settings` (one `name = value` per line, as gfortran writes
it). Later lines win, so a switch row can append overrides before the closing `/`.
"""
function read_settings(path::AbstractString)
    values = Dict{Symbol, Any}()
    for line in eachline(path)
        text = strip(line)
        (isempty(text) || startswith(text, '&') || text == "/") && continue
        name, value = strip.(split(text, '='; limit = 2))
        name = lowercase(name)
        name in NAMELIST_IGNORED && continue
        haskey(NAMELIST_FIELDS, name) || error("unknown namelist entry $name in $path")
        field = NAMELIST_FIELDS[name]
        values[field] = parse_value(fieldtype(Settings, field), rstrip(value, ','))
    end
    return Settings(; values...)
end

function parse_value(::Type{Bool}, text)
    letter = uppercase(strip(text, ['.', ' ']))[1]
    letter in ('T', 'F') || error("not a Fortran logical: $text")
    return letter == 'T'
end
parse_value(::Type{Float64}, text) = parse(Float64, replace(strip(text), r"[dD]" => "e"))
parse_value(::Type{Symbol}, text) = Symbol(parse_value(String, text))
parse_value(::Type{String}, text) = String(strip(strip(text), ['"', '\'', ' ']))
