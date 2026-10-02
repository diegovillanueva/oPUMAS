# 🧊🌧️ Freezing of rain.

"""
    rain_immersion_freezing(parameters, T, q_rain, N_rain, λ_rain; enabled = true)

Heterogeneous (immersion) freezing of raindrops below -4 °C, after Bigg (1953):
the number rate is π N B (e^(A (T₀ - T)) - 1) / λ³, and each frozen drop carries
20 π ρ_w / λ³ (in-precipitation, ≥ 0). Frozen rain becomes snow.
🧊🎛️ Sppe: `enabled = false` (`rain_freeze`) switches it off.
Fortran twin: `heterogeneous_rain_freezing` (`mnuccr`, `nnuccr`).
"""
function rain_immersion_freezing(parameters, T, q_rain, N_rain, λ_rain; enabled = true)
    if !(T < 269.15 && q_rain >= smallest_mixing_ratio && enabled)
        return (; mass = 0.0, number = 0.0)
    end
    (; a, b) = bigg_freezing
    T₀ = parameters.host.melting_point
    number = π * N_rain * b * (exp(a * (T₀ - T)) - 1.0) / λ_rain^3
    mass = number * 20.0 * π * bulk_density.liquid_water / λ_rain^3
    return (; mass, number)
end
