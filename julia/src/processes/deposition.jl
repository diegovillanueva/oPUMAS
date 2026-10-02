# 🌫️❄️ Vapor diffusion: deposition, sublimation, evaporation and the Bergeron process.
#    `air` is a NamedTuple of the ambient air at one point:
#      T (K), ρ (kg m⁻³), diffusivity of vapor (m² s⁻¹), viscosity (kg m⁻¹ s⁻¹),
#      schmidt number (1), q_vapor, q_sat_liquid, q_sat_ice (kg kg⁻¹).

"""
    ventilated_diffusion(air, N₀, λ, a, Γ, b, f)

Vapor exchange coefficient of a gamma population with fall speed V = a Dᵇ (s⁻¹):

    2π N₀ ρ D_v (f₁ / λ² + f₂ (a ρ / μ)^½ Sc^⅓ Γ / λ^(5/2 + b/2))

with ventilation coefficients `f` = (f1, f2) and Γ = Γ(5/2 + b/2). The rate is this
times the supersaturation divided by [`latent_heat_factor`](@ref).
Fortran twin: the `eps` expression repeated in the routines below.
"""
function ventilated_diffusion(air, N₀, λ, a, Γ, b, f)
    (; ρ, diffusivity, viscosity, schmidt) = air
    ventilation_term = f.f2 * (a * ρ / viscosity)^0.5 * schmidt^(1.0 / 3.0) * Γ /
                       λ^(5.0 / 2.0 + b / 2.0)
    return 2.0 * π * N₀ * ρ * diffusivity * (f.f1 / (λ * λ) + ventilation_term)
end

"""
    bergeron_to_snow(parameters, air, a_snow, q_cloud, q_snow, λ_snow, N₀_snow)

Wegener-Bergeron-Findeisen growth of snow at the expense of cloud droplets below 0 °C:
snow grows by deposition at the vapor excess q_sat_liquid - q_sat_ice (kg kg⁻¹ s⁻¹).
Fortran twin: `bergeron_process_snow` (`bergs`).
"""
function bergeron_to_snow(parameters, air, a_snow, q_cloud, q_snow, λ_snow, N₀_snow)
    (; T, q_sat_liquid, q_sat_ice) = air
    if q_snow < smallest_mixing_ratio || q_cloud < smallest_mixing_ratio ||
       T >= parameters.host.melting_point
        return 0.0
    end
    heat = latent_heat_factor(parameters, T, q_sat_ice, parameters.latent_heat_sublimation)
    b = fall_speed_law.snow.b
    exchange = ventilated_diffusion(
        air, N₀_snow, λ_snow, a_snow, parameters.gamma.half_snow_plus5, b, ventilation.snow
    )
    return exchange * (q_sat_liquid - q_sat_ice) / heat
end

"""
    snow_deposition(parameters, air, q_snow, N_snow, precipitation_fraction, a_snow)

Vapor deposition onto snow below 0 °C, for grid-mean snow `q_snow`, `N_snow` spread
over `precipitation_fraction` (only where it is at least 0.1). Grid-mean rate
(kg kg⁻¹ s⁻¹), never negative: snow sublimation is in [`precipitation_evaporation`](@ref).
Fortran twin: `vapor_deposition_onto_snow` (`vap_deps`).
"""
function snow_deposition(parameters, air, q_snow, N_snow, precipitation_fraction, a_snow)
    (; T, q_vapor, q_sat_ice) = air
    if q_snow < smallest_mixing_ratio || precipitation_fraction < 0.1
        return 0.0
    end
    q_in = q_snow / precipitation_fraction
    N_in = N_snow / precipitation_fraction
    heat = latent_heat_factor(parameters, T, q_sat_ice, parameters.latent_heat_sublimation)
    snow = size_distribution(parameters.hydrometeors.snow, q_in, N_in)
    exchange = ventilated_diffusion(air, snow.intercept, snow.slope, a_snow,
        parameters.gamma.half_snow_plus5, fall_speed_law.snow.b, ventilation.snow)
    rate = exchange / heat * (q_vapor - q_sat_ice) * precipitation_fraction
    return T < parameters.host.melting_point && rate > 0.0 ? rate : 0.0
end

"""
    ice_deposition(parameters, air, q_ice, N_ice, ice_fraction)

Vapor growth and sublimation of cloud ice, without ventilation (small crystals), for
grid-mean `q_ice`, `N_ice` spread over `ice_fraction`. Returns grid-mean rates
(kg kg⁻¹ s⁻¹):

- `deposition` ≥ 0, below 0 °C only
- `sublimation` ≤ 0, otherwise
- `bergeron` ≥ 0: in-cloud growth at the vapor excess over ice of liquid saturation,
  used for the Bergeron process on cloud ice

Fortran twin: `ice_deposition_sublimation` (`vap_dep`, `ice_sublim`, `berg`).
"""
function ice_deposition(parameters, air, q_ice, N_ice, ice_fraction)
    (; T, ρ, diffusivity, q_vapor, q_sat_liquid, q_sat_ice) = air
    q_ice >= smallest_mixing_ratio ||
        return (; deposition = 0.0, sublimation = 0.0, bergeron = 0.0)
    heat = latent_heat_factor(parameters, T, q_sat_ice, parameters.latent_heat_sublimation)
    ice = size_distribution(
        parameters.hydrometeors.cloud_ice, q_ice / ice_fraction, N_ice / ice_fraction
    )
    exchange = 2.0 * π * ice.intercept * ρ * diffusivity / (ice.slope * ice.slope)
    growth = exchange / heat * (q_vapor - q_sat_ice) * ice_fraction
    freezing = T < parameters.host.melting_point
    deposition, sublimation = freezing && growth > 0.0 ? (growth, 0.0) :
                              (0.0, min(growth, 0.0))
    bergeron = freezing ? max(exchange / heat * (q_sat_liquid - q_sat_ice), 0.0) : 0.0
    return (; deposition, sublimation, bergeron)
end

"""
    precipitation_evaporation(parameters, air, liquid_fraction, precipitation_fraction,
                              q_cloud, q_ice, rain, snow, graupel)

Evaporation of rain and sublimation of snow and graupel in the clear-air part of the
precipitation area. `rain`, `snow`, `graupel` are NamedTuples `(; q, λ, N₀, a)` of
in-precipitation values (graupel also `b`). Precipitation under cloud (`liquid_fraction`,
ignored if `q_cloud + q_ice < 10⁻⁶`) does not evaporate; vapor in the clear part is
`q_clear = (q_vapor - f_cloud q_sat_liquid) / (1 - f_cloud)`.

Returns in-precipitation rates ≤ 0 (kg kg⁻¹ s⁻¹) `rain`, `snow`, `graupel`, and
`area`, the precipitation fraction outside cloud (`am_evp_st`).
With `ifs_evaporation_threshold`, rain only evaporates where q_clear / q_sat_liquid ≤ 0.9.
Fortran twin: `evaporate_sublimate_precip_graupel` (`pre`, `prds`, `prdg`).
"""
function precipitation_evaporation(
        parameters, air, liquid_fraction, precipitation_fraction,
        q_cloud, q_ice, rain, snow, graupel)
    (; T, q_vapor, q_sat_liquid, q_sat_ice) = air
    none = (; rain = 0.0, snow = 0.0, graupel = 0.0, area = 0.0)
    cloud_fraction = q_cloud + q_ice < 1.0e-6 ? 0.0 : liquid_fraction
    precipitation_fraction > cloud_fraction || return none
    small = smallest_mixing_ratio
    rain.q >= small || snow.q >= small || graupel.q >= small || return none
    area = precipitation_fraction - cloud_fraction
    q_clear = (q_vapor - cloud_fraction * q_sat_liquid) / (1.0 - cloud_fraction)
    # 🧮 grid-mean loss, capped at zero, back to in-precipitation
    in_precipitation(rate) = min(rate * area, 0.0) / precipitation_fraction
    L_s = parameters.latent_heat_sublimation
    b_snow = fall_speed_law.snow.b

    rain_rate = 0.0
    dry_enough = !parameters.settings.ifs_evaporation_threshold ||
                 q_clear / q_sat_liquid <= 0.9
    if rain.q >= small && dry_enough
        L_v = parameters.host.latent_heat_vaporization
        heat = latent_heat_factor(parameters, T, q_sat_liquid, L_v)
        exchange = ventilated_diffusion(air, rain.N₀, rain.λ, rain.a,
            parameters.gamma.half_rain_plus5, fall_speed_law.rain.b, ventilation.rain)
        rain_rate = in_precipitation(exchange * (q_clear - q_sat_liquid) / heat)
    end

    snow_rate = 0.0
    if snow.q >= small
        heat = latent_heat_factor(parameters, T, q_sat_ice, L_s)
        exchange = ventilated_diffusion(air, snow.N₀, snow.λ, snow.a,
            parameters.gamma.half_snow_plus5, b_snow, ventilation.snow)
        snow_rate = in_precipitation(exchange * (q_clear - q_sat_ice) / heat)
    end

    graupel_rate = 0.0
    if graupel.q >= small
        heat = latent_heat_factor(parameters, T, q_sat_ice, L_s)
        # 🐛 Fortran uses Γ(5/2 + b_g/2) but λ^(5/2 + b_s/2): see FINDINGS.md.
        Γ = gamma_function(5.0 / 2.0 + graupel.b / 2.0)
        exchange = ventilated_diffusion(
            air, graupel.N₀, graupel.λ, graupel.a, Γ, b_snow, ventilation.snow
        )
        graupel_rate = in_precipitation(exchange * (q_clear - q_sat_ice) / heat)
    end
    return (; rain = rain_rate, snow = snow_rate, graupel = graupel_rate, area)
end
