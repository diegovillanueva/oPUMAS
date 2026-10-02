# 🟤🧊 Graupel: collection by graupel, conversion of rimed snow, rime splintering.
#    Collision kernels of Ikawa and Saito (1991) and Reisner et al. (1998).

"""
    snow_collection_by_rain(U_rain_mass, U_snow_mass, ρ_air, q_snow, q_rain,
                            λ_rain, N₀_rain, λ_snow, N₀_snow)

Snow collected by rain, which becomes graupel (in-precipitation, kg kg⁻¹ s⁻¹, ≥ 0).
Applies only when both snow and rain exceed 0.1 g kg⁻¹.
Fortran twin: `graupel_collecting_snow` (`psacr`).
"""
function snow_collection_by_rain(U_rain_mass, U_snow_mass, ρ_air, q_snow, q_rain,
        λ_rain, N₀_rain, λ_snow, N₀_snow)
    q_snow >= 0.1e-3 && q_rain >= 0.1e-3 || return 0.0
    kernel = π * π * collection_efficiency.droplet_rain * bulk_density.snow
    speed = ((1.2 * U_rain_mass - 0.95 * U_snow_mass)^2 +
             0.08 * U_snow_mass * U_rain_mass)^0.5
    moments = 5.0 / (λ_snow^3 * λ_rain) + 2.0 / (λ_snow^2 * λ_rain^2) +
              0.5 / (λ_snow * λ_rain^3)
    return kernel * (speed * ρ_air * N₀_rain * N₀_snow / λ_snow^3 * moments)
end

"""
    droplet_collection_by_graupel(q_graupel, q_cloud, N_cloud, ρ_air, N₀_graupel,
                                  λ_graupel, b_graupel, a_graupel)

Riming of graupel: collection of cloud droplets by graupel with fall speed a Dᵇ
(in-cloud, ≥ 0). Applies when graupel exceeds `small_incloud_mixing_ratio`.
Fortran twin: `graupel_collecting_cld_water` (`psacwg`, `npsacwg`).
"""
function droplet_collection_by_graupel(q_graupel, q_cloud, N_cloud, ρ_air, N₀_graupel,
        λ_graupel, b_graupel, a_graupel)
    if !(q_graupel >= 1.0e-8 && q_cloud >= smallest_mixing_ratio)
        return (; mass = 0.0, number = 0.0)
    end
    kernel = gamma_function(b_graupel + 3.0) * π / 4.0 * collection_efficiency.ice_droplet
    λ_power = λ_graupel^(b_graupel + 3.0)
    mass = kernel * a_graupel * q_cloud * ρ_air * N₀_graupel / λ_power
    number = kernel * a_graupel * N_cloud * ρ_air * N₀_graupel / λ_power
    return (; mass, number)
end

"""
    rimed_snow_to_graupel(parameters, riming, q_snow, q_cloud, N_snow, ρ_air,
                          snow_density, graupel_density, a_snow, λ_snow, N₀_snow,
                          time_step)

Heavily rimed snow becomes graupel (Reisner et al. 1998): part of the snow riming
`riming` (kg kg⁻¹ s⁻¹) goes to graupel when riming > 0, q_snow ≥ 0.1 g kg⁻¹ and
q_cloud ≥ 0.5 g kg⁻¹. Returns the graupel `mass` rate, the snow `number` converted
to graupel particles of `new_graupel_mass` (at most N_snow per time step), and the
`riming` left for snow.
Fortran twin: `graupel_riming_liquid_snow` (`pgsacw`, `nscng`, `psacws`).
"""
function rimed_snow_to_graupel(parameters, riming, q_snow, q_cloud, N_snow, ρ_air,
        snow_density, graupel_density, a_snow, λ_snow, N₀_snow, time_step)
    if !(riming > 0.0 && q_snow >= 0.1e-3 && q_cloud >= 0.5e-3)
        return (; mass = 0.0, number = 0.0, riming)
    end
    E = collection_efficiency.ice_droplet
    b = fall_speed_law.snow.b
    ρ_reference = parameters.reference_air_density
    factor = 4.0 * 2.0 * 3.0 * ρ_reference * π * E * E *
             parameters.gamma.twice_snow_plus2 / (8.0 * (graupel_density - snow_density))
    mass = min(riming,
        factor * time_step * N₀_snow * q_cloud * q_cloud *
        a_snow * a_snow / (ρ_air * λ_snow^(2.0 * b + 2.0)))
    volume_gain = max(snow_density / (graupel_density - snow_density) * mass, 0.0)
    number = volume_gain / new_graupel_mass * ρ_air
    number = min(number, N_snow / time_step)
    return (; mass, number, riming = riming - mass)
end

"""
    rain_collection_by_graupel(q_rain, q_graupel, U_graupel_mass, U_rain_mass,
                               U_graupel_number, U_rain_number, ρ_air,
                               N₀_rain, λ_rain, N₀_graupel, λ_graupel)

Collection of rain by graupel (in-precipitation, ≥ 0 for mass). The number rate
removes the collected drops and adds back the mass at 5.2 × 10⁻⁷ kg per drop.
Applies when rain and graupel exceed `small_incloud_mixing_ratio`.
Fortran twin: `graupel_collecting_rain` (`pracg`, `npracg`).
"""
function rain_collection_by_graupel(q_rain, q_graupel, U_graupel_mass, U_rain_mass,
        U_graupel_number, U_rain_number, ρ_air, N₀_rain, λ_rain, N₀_graupel, λ_graupel)
    q_rain >= 1.0e-8 && q_graupel >= 1.0e-8 || return (; mass = 0.0, number = 0.0)
    E = collection_efficiency.droplet_rain
    mass_kernel = π * π * E * bulk_density.liquid_water
    number_kernel = π / 2.0 * E
    # 🧮 The Fortran writes real exponents (**3._r8) here, unlike in the snow kernels.
    λr, λg = λ_rain, λ_graupel
    mass_speed = ((1.2 * U_rain_mass - 0.95 * U_graupel_mass)^2.0 +
                  0.08 * U_graupel_mass * U_rain_mass)^0.5
    mass_moments = 5.0 / (λr^3.0 * λg) + 2.0 / (λr^2.0 * λg^2.0) + 0.5 / (λr * λg^3.0)
    mass = mass_kernel *
           (mass_speed * ρ_air * N₀_rain * N₀_graupel / λr^3.0 * mass_moments)
    number_speed = (1.7 * (U_rain_number - U_graupel_number)^2.0 +
                    0.3 * U_rain_number * U_graupel_number)^0.5
    number_moments = 1.0 / (λr^3.0 * λg) + 1.0 / (λr^2.0 * λg^2.0) + 1.0 / (λr * λg^3.0)
    number = number_kernel * ρ_air * number_speed * N₀_rain * N₀_graupel * number_moments
    return (; mass, number = number - mass / 5.2e-7)
end

"""
    rain_snow_collisions_to_graupel(rain_by_snow, rain_number_by_snow, snow_by_rain,
                                    q_snow, q_rain, N_rain, N_snow, N₀_snow, λ_snow,
                                    N₀_rain, λ_rain, time_step)

Rain-snow collisions make graupel (Reisner et al. 1998). The fraction that stays
snow is f = ρ_s² (4/λ_s)⁶ / (ρ_s² (4/λ_s)⁶ + ρ_w² (4/λ_r)⁶), the mass-weighted share
of snow. The rest of rain collected by snow (`rain_by_snow`, number
`rain_number_by_snow`) goes to graupel, as does the rest of `snow_by_rain`.
Applies when rain_by_snow > 0 and both snow and rain exceed 0.1 g kg⁻¹.

Returns the graupel `mass` and `number` rates and what is left of the three inputs.
Fortran twin: `graupel_rain_riming_snow` (`pgracs`, `ngracs`, `pracs`, `npracs`, `psacr`).
"""
function rain_snow_collisions_to_graupel(rain_by_snow, rain_number_by_snow, snow_by_rain,
        q_snow, q_rain, N_rain, N_snow, N₀_snow, λ_snow, N₀_rain, λ_rain, time_step)
    if !(rain_by_snow > 0.0 && q_snow >= 0.1e-3 && q_rain >= 0.1e-3)
        return (; mass = 0.0, number = 0.0, rain_by_snow, rain_number_by_snow,
            snow_by_rain)
    end
    ρ_s², ρ_w² = bulk_density.snow^2, bulk_density.liquid_water^2
    snow_weight = ρ_s² * (4.0 / λ_snow)^3 * (4.0 / λ_snow)^3
    rain_weight = ρ_w² * (4.0 / λ_rain)^3 * (4.0 / λ_rain)^3
    stays_snow = snow_weight / (ρ_s² * (4.0 / λ_snow)^3 * (4.0 / λ_snow)^3 + rain_weight)
    stays_snow = max(min(stays_snow, 1.0), 0.0)
    mass = (1.0 - stays_snow) * rain_by_snow
    number = (1.0 - stays_snow) * rain_number_by_snow
    number = min(number, N_rain / time_step)
    number = min(number, N_snow / time_step)
    return (; mass, number,
        rain_by_snow = rain_by_snow - mass,
        rain_number_by_snow = rain_number_by_snow - number,
        snow_by_rain = snow_by_rain * (1.0 - stays_snow))
end

"""
    rime_splintering_by_graupel(parameters, T, q_cloud, q_rain, q_graupel,
                                droplet_riming, rain_riming)

Hallett-Mossop rime splintering on graupel (graupel ≥ 0.1 g kg⁻¹ riming cloud
water ≥ 0.5 g kg⁻¹ or rain ≥ 0.1 g kg⁻¹) between -3 and -8 °C: 3.5 × 10⁸ splinters
per kg of rime, peaking at -5 °C. Splinters of `splinter_mass` are taken from the
rime of droplets (`droplet_riming`) and of rain (`rain_riming`).

Returns the ice `number` and `mass` from each source and the rime rates left.
Fortran twin: `graupel_rime_splintering` (`nmultg`, `qmultg`, `nmultrg`, `qmultrg`).
"""
function rime_splintering_by_graupel(parameters, T, q_cloud, q_rain, q_graupel,
        droplet_riming, rain_riming)
    T₀ = parameters.host.melting_point
    T_3, T_5, T_8 = T₀ - 3.0, T₀ - 5.0, T₀ - 8.0
    none = (; droplet_mass = 0.0, droplet_number = 0.0, rain_mass = 0.0,
        rain_number = 0.0, droplet_riming, rain_riming)
    q_graupel >= 0.1e-3 || return none
    q_cloud >= 0.5e-3 || q_rain >= 0.1e-3 || return none
    droplet_riming > 0.0 || rain_riming > 0.0 || return none
    T < T_3 && T > T_8 || return none
    # 🐛 The Fortran also tests T > T_3 and T < T_8 here, which cannot happen.
    efficiency = T > T_5 ? (T_3 - T) / 2.0 : (T - T_8) / 3.0

    droplet_mass, droplet_number = 0.0, 0.0
    if droplet_riming > 0.0
        droplet_number = 35.0e4 * droplet_riming * efficiency * 1000.0
        droplet_mass = min(droplet_number * splinter_mass, droplet_riming)
        droplet_riming = droplet_riming - droplet_mass
    end
    rain_mass, rain_number = 0.0, 0.0
    if rain_riming > 0.0
        rain_number = 35.0e4 * rain_riming * efficiency * 1000.0
        rain_mass = min(rain_number * splinter_mass, rain_riming)
        rain_riming = rain_riming - rain_mass
    end
    return (; droplet_mass, droplet_number, rain_mass, rain_number, droplet_riming,
        rain_riming)
end
