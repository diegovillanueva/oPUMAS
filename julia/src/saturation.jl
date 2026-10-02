# 💧🌡️ Saturation vapor pressure (Goff and Gratch, 1946) and saturation specific humidity.

"""
Constants of the saturation routines, as given to `wv_sat_methods_init`.

- `melting_point`: K (`tmelt`)
- `triple_point`: K (`h2otrip`)
- `boiling_point`: steam point, K (`tboil`)
- `ice_transition_range`: K (`ttrice`, unused by Goff-Gratch)
- `molar_mass_ratio`: water vapor over dry air, ε (`epsilo`)
"""
Base.@kwdef struct SaturationConstants
    melting_point::Float64
    triple_point::Float64
    boiling_point::Float64
    ice_transition_range::Float64
    molar_mass_ratio::Float64
end

"""
    saturation_vapor_pressure_over_water(sat, T)

Goff-Gratch saturation vapor pressure over liquid water (Pa) at temperature `T` (K).
Uncertain below -70 °C.

Fortran twin: `GoffGratch_svp_water` in `wv_sat_methods.F90`.
"""
function saturation_vapor_pressure_over_water(sat::SaturationConstants, T)
    Tₛ = sat.boiling_point
    log10_e_hPa = -7.90298 * (Tₛ / T - 1.0) +
                  5.02808 * log10(Tₛ / T) -
                  1.3816e-7 * (10.0^(11.344 * (1.0 - T / Tₛ)) - 1.0) +
                  8.1328e-3 * (10.0^(-3.49149 * (Tₛ / T - 1.0)) - 1.0) +
                  log10(1013.246)
    return 10.0^log10_e_hPa * 100.0
end

"""
    saturation_vapor_pressure_over_ice(sat, T)

Goff-Gratch saturation vapor pressure over ice (Pa) at temperature `T` (K).
Good down to -100 °C.

Fortran twin: `GoffGratch_svp_ice` in `wv_sat_methods.F90`.
"""
function saturation_vapor_pressure_over_ice(sat::SaturationConstants, T)
    T₀ = sat.triple_point
    log10_e_hPa = -9.09718 * (T₀ / T - 1.0) - 3.56654 * log10(T₀ / T) +
                  0.876793 * (1.0 - T / T₀) +
                  log10(6.1071)
    return 10.0^log10_e_hPa * 100.0
end

"""
    vapor_pressure_to_specific_humidity(sat, e, p)

Saturation specific humidity (kg kg⁻¹) from saturation vapor pressure `e` and air
pressure `p` (both Pa). Capped at 1 where `p ≤ e`.

Fortran twin: `wv_sat_svp_to_qsat` in `wv_sat_methods.F90`.
"""
function vapor_pressure_to_specific_humidity(sat::SaturationConstants, e, p)
    ε = sat.molar_mass_ratio
    p - e <= 0.0 && return 1.0
    return ε * e / (p - (1.0 - ε) * e)
end

"""
    saturation_over_water(sat, T, p) -> (; vapor_pressure, specific_humidity)

Saturation vapor pressure (Pa, capped at `p`) and specific humidity (kg kg⁻¹) over
liquid water. Fortran twin: `wv_sat_qsat_water` (as `qsat_water`) in `wv_sat_methods.F90`.
"""
function saturation_over_water(sat::SaturationConstants, T, p)
    e = saturation_vapor_pressure_over_water(sat, T)
    q = vapor_pressure_to_specific_humidity(sat, e, p)
    return (; vapor_pressure = min(e, p), specific_humidity = q)
end

"""
    saturation_over_ice(sat, T, p) -> (; vapor_pressure, specific_humidity)

As [`saturation_over_water`](@ref), over ice.
Fortran twin: `wv_sat_qsat_ice` (as `qsat_ice`) in `wv_sat_methods.F90`.
"""
function saturation_over_ice(sat::SaturationConstants, T, p)
    e = saturation_vapor_pressure_over_ice(sat, T)
    q = vapor_pressure_to_specific_humidity(sat, e, p)
    return (; vapor_pressure = min(e, p), specific_humidity = q)
end
