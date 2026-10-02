# ❄️🤝 Ice and snow: autoconversion of cloud ice to snow, aggregation, collection of ice.

"""
    ice_to_snow_autoconversion(parameters, T, q_ice, λ_ice, N₀_ice)

Cloud ice crystals larger than the autoconversion size D_cs become snow within
τ = 180 s, at or below 0 °C. With x = λ D_cs, the number rate is
N₀ / (λ τ) e^(-x) and the mass rate is that times (ρ_i π / 6) / λ³ · (x³ + 3x² + 6x + 6)
(in-cloud, kg kg⁻¹ s⁻¹ and kg⁻¹ s⁻¹, both ≥ 0).
Fortran twin: `ice_autoconversion` (`prci`, `nprci`).
"""
function ice_to_snow_autoconversion(parameters, T, q_ice, λ_ice, N₀_ice)
    if !(T <= parameters.host.melting_point && q_ice >= smallest_mixing_ratio)
        return (; mass = 0.0, number = 0.0)
    end
    conversion_time = 180.0
    x = λ_ice * parameters.settings.ice_to_snow_size
    number = N₀_ice / (λ_ice * conversion_time) * exp(-x)
    mean_mass = (bulk_density.cloud_ice * π / 6.0) / λ_ice^3
    mass = mean_mass * number * (((x + 3.0) * x + 6.0) * x + 6.0)
    return (; mass, number)
end

"""
    snow_self_aggregation(parameters, T, ρ_air, a_snow, snow_density, q_snow, N_snow)

Loss of snow number by snow flakes sticking together, at or below 0 °C (kg⁻¹ s⁻¹, ≤ 0).
Fortran twin: `snow_self_aggregation` (`nsagg`).
"""
function snow_self_aggregation(parameters, T, ρ_air, a_snow, snow_density, q_snow,
        N_snow)
    if !(q_snow >= smallest_mixing_ratio && T <= parameters.host.melting_point)
        return 0.0
    end
    E = collection_efficiency.ice_snow
    b = fall_speed_law.snow.b
    mean_mass_per_volume = (q_snow / N_snow) * (1.0 / (snow_density * π))
    return -(1108.0 * E / (4.0 * 720.0 * snow_density) * a_snow * q_snow * N_snow *
             ρ_air * mean_mass_per_volume^((b - 1.0) / 3.0))
end

"""
    ice_collection_by_snow(parameters, T, ρ_air, a_snow, q_ice, N_ice, q_snow,
                           λ_snow, N₀_snow)

Collection of cloud ice by falling snow at or below 0 °C (in-cloud, ≥ 0). The kernel
π/4 E a ρ N₀ Γ(3+b) / λ^(3+b) is applied to both mass and number of the ice.
Fortran twin: `accrete_cloud_ice_snow` (`prai`, `nprai`).
"""
function ice_collection_by_snow(parameters, T, ρ_air, a_snow, q_ice, N_ice, q_snow,
        λ_snow, N₀_snow)
    small = smallest_mixing_ratio
    if !(q_snow >= small && q_ice >= small && T <= parameters.host.melting_point)
        return (; mass = 0.0, number = 0.0)
    end
    b = fall_speed_law.snow.b
    kernel = π / 4.0 * collection_efficiency.ice_snow * a_snow * ρ_air * N₀_snow *
             parameters.gamma.snow_plus3 / λ_snow^(b + 3.0)
    return (; mass = kernel * q_ice, number = kernel * N_ice)
end
