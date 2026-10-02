# 🎛️🧭 Settings (the namelist) and Parameters (everything micro_pumas_init computes).

"""
Physical constants of the host model, passed to `micro_pumas_init`. Defaults are the
CESM `shr_const` values.

- `gravity` (m s⁻², `gravit`)
- `gas_constant_dry_air`, `gas_constant_water_vapor` (J kg⁻¹ K⁻¹, `rair`, `rh2o`)
- `specific_heat_dry_air` (J kg⁻¹ K⁻¹, `cpair`)
- `melting_point`, `triple_point`, `boiling_point` (K, `tmelt`, `h2otrip`, `tboil`)
- `latent_heat_vaporization`, `latent_heat_fusion` (J kg⁻¹, `latvap`, `latice`)
- `ice_transition_range` (K, `ttrice`)
- `molar_mass_ratio` (water vapor over dry air, `epsilo`)
"""
Base.@kwdef struct HostConstants
    gravity::Float64 = 9.80616
    gas_constant_dry_air::Float64 = 6.02214e26 * 1.38065e-23 / 28.966
    gas_constant_water_vapor::Float64 = 6.02214e26 * 1.38065e-23 / 18.016
    specific_heat_dry_air::Float64 = 1.00464e3
    melting_point::Float64 = 273.15
    triple_point::Float64 = 273.16
    boiling_point::Float64 = 373.16
    latent_heat_vaporization::Float64 = 2.501e6
    latent_heat_fusion::Float64 = 3.337e5
    ice_transition_range::Float64 = 20.0
    molar_mass_ratio::Float64 = 18.016 / 28.966
end

"""
The PUMAS namelist. Each field names its CAM namelist variable. Defaults are the
upstream CESM 2.3 behaviour (CESM 2.1.5 interface defaults, CAM6 with CLUBB).
"""
Base.@kwdef struct Settings
    "Smallest relative humidity for ice cloud fraction > 0 (`cldfrc2m_rhmini`)."
    minimum_rh_for_ice_cloud::Float64 = 0.80
    "Cloud ice to snow autoconversion diameter, m (`micro_mg_dcs`)."
    ice_to_snow_size::Float64 = 500.0e-6
    "Prognostic hail instead of graupel (`micro_mg_do_hail`)."
    do_hail::Bool = false
    "Prognostic graupel (`micro_mg_do_graupel`)."
    do_graupel::Bool = false
    "Sub-columns are uniform (`microp_uniform`)."
    uniform_subcolumns::Bool = false
    "Cloud ice is prognostic (`do_cldice`)."
    prognostic_cloud_ice::Bool = true
    "Classical-nucleation heterogeneous freezing from the host (`use_hetfrz_classnuc`)."
    host_heterogeneous_freezing::Bool = true
    "Precipitation fraction: `:in_cloud` or `:max_overlap` (`micro_mg_precip_frac_method`)."
    precipitation_fraction_method::Symbol = :in_cloud
    "Bergeron efficiency factor (`micro_mg_berg_eff_factor`)."
    bergeron_efficiency::Float64 = 1.0
    "Accretion enhancement factor (`micro_mg_accre_enhan_fact`)."
    accretion_enhancement::Float64 = 1.0
    "Autoconversion prefactor (`micro_mg_autocon_fact`)."
    autoconversion_prefactor::Float64 = 0.01
    "Autoconversion exponent of droplet number (`micro_mg_autocon_nd_exp`)."
    autoconversion_number_exponent::Float64 = -1.1
    "Autoconversion exponent of liquid water (`micro_mg_autocon_lwp_exp`)."
    autoconversion_water_exponent::Float64 = 2.47
    "Diameter of homogeneously frozen ice, m (`micro_mg_homog_size`)."
    homogeneous_freezing_size::Float64 = 25.0e-6
    "Factor on the cloud ice fall speed (`micro_mg_vtrmi_factor`)."
    ice_fall_speed_factor::Float64 = 1.0
    "Factor on the snow fall speed (`micro_mg_vtrms_factor`)."
    snow_fall_speed_factor::Float64 = 1.0
    "Factor on the ice effective radius (`micro_mg_effi_factor`)."
    ice_effective_radius_factor::Float64 = 1.0
    "Factor on accretion of cloud ice by snow (`micro_mg_iaccr_factor`)."
    ice_accretion_factor::Float64 = 1.0
    "Largest ice number concentration, m⁻³ (`micro_mg_max_nicons`)."
    maximum_ice_number::Float64 = 1.0e8
    """
    Remove supersaturation after sedimentation (`remove_supersat` in PUMAS). CAM passes
    its namelist variable `allow_sed_supersat` here unchanged: see FINDINGS.md.
    """
    remove_supersaturation::Bool = false
    "Warm rain scheme: `:kk2000` or `:sb2001` (`micro_pumas_warm_rain`)."
    warm_rain::Symbol = :kk2000
    "No evaporation of sedimenting condensate by cloud fraction (`micro_mg_evap_sed_off`)."
    no_sedimentation_evaporation::Bool = false
    "Ice nucleation without the humidity condition (`micro_mg_icenuc_rh_off`)."
    ice_nucleation_without_rh::Bool = false
    "Meyers et al. (1992) ice nucleation (`micro_mg_icenuc_use_meyers`)."
    meyers_ice_nucleation::Bool = false
    "Scale evaporation by 0.3 as in IFS (`micro_mg_evap_scl_ifs`)."
    ifs_evaporation_scaling::Bool = false
    "Evaporation humidity threshold as in IFS (`micro_mg_evap_rhthrsh_ifs`)."
    ifs_evaporation_threshold::Bool = false
    "Rain freezes at 0 °C as in IFS (`micro_mg_rainfreeze_ifs`)."
    ifs_rain_freezing::Bool = false
    "Snow falls at 1 m s⁻¹ as in IFS (`micro_mg_ifs_sed`)."
    ifs_snow_fall_speed::Bool = false
    "Rain fall speed stays non-zero below rain (`micro_mg_precip_fall_corr`)."
    precipitation_fall_correction::Bool = false
    "Accretion sees the rain made by autoconversion (`micro_mg_accre_sees_auto`)."
    accretion_sees_autoconversion::Bool = false
    "Implicit sedimentation (`micro_mg_implicit_fall`)."
    implicit_sedimentation::Bool = false
    "Hold droplet number constant (`micro_mg_nccons`)."
    constant_droplet_number::Bool = false
    "Hold ice number constant (`micro_mg_nicons`)."
    constant_ice_number::Bool = false
    "Droplet number when held constant, m⁻³ (`micro_mg_ncnst`)."
    droplet_number::Float64 = 100.0e6
    "Ice number when held constant, m⁻³ (`micro_mg_ninst`)."
    ice_number::Float64 = 0.1e6
    "Hold graupel number constant (`micro_mg_ngcons`)."
    constant_graupel_number::Bool = false
    "Graupel number when held constant, m⁻³ (`micro_mg_ngnst`)."
    graupel_number::Float64 = 0.0005e6
end

"""
Everything `micro_pumas_init` and `micro_pumas_utils_init` store, built by [`setup`](@ref).

- `settings`, `host`: the inputs
- `saturation`: constants of the saturation routines
- `hydrometeors`: [`HydrometeorProperties`](@ref) of the six categories
- `latent_heat_sublimation`: L_v + L_f (`xxls`)
- `reference_air_density`: air density at 850 hPa and 0 °C, kg m⁻³ (`rhosu`)
- `ice_nucleation_temperature`, `snow_melting_temperature`, `rain_freezing_temperature`:
  K (`icenuct`, `snowmelt`, `rainfrze`)
- `gamma`: Γ values of the fall speed exponents (`gamma_b*_plus*`, `gamma_*bs*`)
"""
struct Parameters
    settings::Settings
    host::HostConstants
    saturation::SaturationConstants
    hydrometeors::NamedTuple
    latent_heat_sublimation::Float64
    reference_air_density::Float64
    ice_nucleation_temperature::Float64
    snow_melting_temperature::Float64
    rain_freezing_temperature::Float64
    gamma::NamedTuple
end

"""
    setup(; host = HostConstants(), settings...) -> Parameters

Build the parameters from host constants and namelist settings (keywords of
[`Settings`](@ref)). Fortran twin: `micro_pumas_init` in `micro_pumas_v1.F90`, which
calls `micro_pumas_utils_init` in `micro_pumas_utils.F90`.
"""
setup(; host = HostConstants(), settings...) = setup(Settings(; settings...), host)

function setup(settings::Settings, host::HostConstants = HostConstants())
    check_scope(settings)
    T₀ = host.melting_point
    gas_constant = host.gas_constant_dry_air
    bᵣ, bₛ = fall_speed_law.rain.b, fall_speed_law.snow.b
    bᵢ, bⱼ = fall_speed_law.cloud_ice.b, fall_speed_law.small_cloud_ice.b
    b_dense = dense_precipitation_exponent(settings)
    gamma = (
        snow_plus3 = gamma_function(3.0 + bₛ),
        half_rain_plus5 = gamma_function(5.0 / 2.0 + bᵣ / 2.0),
        half_snow_plus5 = gamma_function(5.0 / 2.0 + bₛ / 2.0),
        twice_snow_plus2 = gamma_function(2.0 * bₛ + 2.0),
        rain_plus1 = gamma_function(1.0 + bᵣ),
        rain_plus4 = gamma_function(4.0 + bᵣ),
        snow_plus1 = gamma_function(1.0 + bₛ),
        snow_plus4 = gamma_function(4.0 + bₛ),
        ice_plus1 = gamma_function(1.0 + bᵢ),
        ice_plus4 = gamma_function(4.0 + bᵢ),
        small_ice_plus1 = gamma_function(1.0 + bⱼ),
        small_ice_plus4 = gamma_function(4.0 + bⱼ),
        graupel_plus1 = gamma_function(1.0 + b_dense),
        graupel_plus4 = gamma_function(4.0 + b_dense),
    )
    saturation = SaturationConstants(;
        melting_point = T₀,
        triple_point = host.triple_point,
        boiling_point = host.boiling_point,
        ice_transition_range = host.ice_transition_range,
        molar_mass_ratio = host.molar_mass_ratio,
    )
    rain_freezing_temperature = settings.ifs_rain_freezing ? T₀ : T₀ - 40.0
    return Parameters(
        settings,
        host,
        saturation,
        hydrometeor_properties(settings.ice_to_snow_size),
        host.latent_heat_vaporization + host.latent_heat_fusion,
        85000.0 / (gas_constant * T₀),
        T₀ - 5.0,
        T₀ + 2.0,
        rain_freezing_temperature,
        gamma,
    )
end

"""
Fall speed exponent of the dense precipitation category: graupel, hail, or 0 if neither
is on (the Fortran then computes Γ(1) and Γ(4)).
"""
function dense_precipitation_exponent(settings::Settings)
    settings.do_graupel && return fall_speed_law.graupel.b
    settings.do_hail && return fall_speed_law.hail.b
    return 0.0
end

"Stop with a clear message for settings that the twin does not translate."
function check_scope(settings::Settings)
    if settings.do_hail && settings.do_graupel
        error("do_hail and do_graupel cannot both be true")
    end
    if settings.warm_rain ∉ (:kk2000, :sb2001)
        error("warm_rain = :$(settings.warm_rain) is not translated yet")
    end
    if settings.precipitation_fraction_method ∉ (:in_cloud, :max_overlap)
        method = settings.precipitation_fraction_method
        error("unknown precipitation_fraction_method $method")
    end
    return nothing
end
