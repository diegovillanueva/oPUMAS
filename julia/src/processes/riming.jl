# ❄️💧 Riming: snow collects cloud droplets and rain; rime splintering makes new ice.

"""
    droplet_collection_by_snow(parameters, T, ρ_air, a_snow, U_snow, viscosity,
                               q_cloud, N_cloud, q_snow, μ_cloud, λ_cloud,
                               λ_snow, N₀_snow)

Riming: collection of cloud droplets by snow at or below 0 °C (in-cloud, ≥ 0).
The collection efficiency follows the Stokes number of the mean droplet,
D = (μ + 1) / λ, falling past snow with number-weighted fall speed `U_snow`:
S = D² U ρ_w λ_s / (9 η), E = S² / (S + 0.4)², limited to [0, 1].
Fortran twin: `accrete_cloud_water_snow` (`psacws`, `npsacws`).
"""
function droplet_collection_by_snow(parameters, T, ρ_air, a_snow, U_snow, viscosity,
        q_cloud, N_cloud, q_snow, μ_cloud, λ_cloud, λ_snow, N₀_snow)
    small = smallest_mixing_ratio
    if !(q_snow >= small && T <= parameters.host.melting_point && q_cloud >= small)
        return (; mass = 0.0, number = 0.0)
    end
    D = (μ_cloud + 1.0) / λ_cloud
    stokes = D * D * U_snow * bulk_density.liquid_water * λ_snow / (9.0 * viscosity)
    efficiency = stokes * stokes / ((stokes + 0.4) * (stokes + 0.4))
    efficiency = min(max(efficiency, 0.0), 1.0)
    b = fall_speed_law.snow.b
    kernel = π / 4.0 * a_snow * ρ_air * N₀_snow * efficiency *
             parameters.gamma.snow_plus3 / λ_snow^(b + 3.0)
    return (; mass = kernel * q_cloud, number = kernel * N_cloud)
end

"""
    rime_splintering_by_snow(T, riming; enabled = true, number_factor = 1.0)

Hallett-Mossop rime splintering: riming snow (`riming`, kg kg⁻¹ s⁻¹) sheds
3.5 × 10⁸ ice splinters per kg of rime, most at -5 °C, none outside -3 to -8 °C.
Returns the new ice `number` (kg⁻¹ s⁻¹), its `mass` (splinters of `new_ice_crystal_mass`,
at most the rime), and the `riming` left for snow.
🧊🎛️ Sppe: `enabled = false` (`onoff_sip`) gives no splinters, `number_factor`
(`factor_HM`) scales their number.
Fortran twin: `secondary_ice_production` (`nsacwi`, `msacwi`, `psacws`).
"""
function rime_splintering_by_snow(T, riming; enabled = true, number_factor = 1.0)
    if !enabled
        number = 0.0
    elseif T < 270.16 && T >= 268.16
        number = 3.5e8 * (270.16 - T) / 2.0 * riming
    elseif T < 268.16 && T >= 265.16
        number = 3.5e8 * (T - 265.16) / 3.0 * riming
    else
        number = 0.0
    end
    number_factor != 1.0 && (number = number * number_factor)
    mass = min(number * new_ice_crystal_mass, riming)
    return (; mass, number, riming = riming - mass)
end

"""
    rain_collection_by_snow(parameters, T, ρ_air, U_rain_mass, U_snow_mass,
                            U_rain_number, U_snow_number, q_rain, q_snow,
                            λ_rain, N₀_rain, λ_snow, N₀_snow)

Collection of rain by snow at or below 0 °C (in-precipitation, ≥ 0), with the
fall speed difference of Mizuno (1990) for mass- and number-weighted speeds U.
Applies only when both rain and snow exceed `small_incloud_mixing_ratio`.
Fortran twin: `accrete_rain_snow` (`pracs`, `npracs`).
"""
function rain_collection_by_snow(parameters, T, ρ_air, U_rain_mass, U_snow_mass,
        U_rain_number, U_snow_number, q_rain, q_snow, λ_rain, N₀_rain, λ_snow, N₀_snow)
    small = small_incloud_mixing_ratio
    if !(q_rain >= small && q_snow >= small && T <= parameters.host.melting_point)
        return (; mass = 0.0, number = 0.0)
    end
    E = collection_efficiency.droplet_rain
    common = π * E * ρ_air * N₀_rain * N₀_snow / (λ_rain^3 * λ_snow)
    r = λ_rain / λ_snow
    mass_speed = sqrt((1.2 * U_rain_mass - 0.95 * U_snow_mass)^2 +
                      0.08 * U_snow_mass * U_rain_mass)
    number_speed = sqrt(1.7 * (U_rain_number - U_snow_number)^2 +
                        0.3 * U_rain_number * U_snow_number)
    mass = common * π * bulk_density.liquid_water * mass_speed / λ_rain^3 *
           ((0.5 * r + 2.0) * r + 5.0)
    number = common * 0.5 * number_speed * ((r + 1.0) * r + 1.0)
    return (; mass, number)
end
