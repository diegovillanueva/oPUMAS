# 🧊💧 Properties of each hydrometeor category used by the size distributions.

"""
Fixed properties of one hydrometeor category (Fortran type `MGHydrometeorProps`).

- `density`: bulk density ρ (kg m⁻³)
- `mass_exponent`: exponent d of the mass-diameter law m = c Dᵈ (`eff_dim`)
- `mass_coefficient`: c = ρ π Γ(d + 1) / 6 (`shape_coef`)
- `slope_bounds`: (smallest, largest) slope λ (m⁻¹), or `nothing` if λ is not limited
- `minimum_mean_mass`: smallest mean particle mass (kg), or `nothing` if not limited

Fortran twin: `NewMGHydrometeorProps` in `micro_pumas_utils.F90`. The Fortran marks a
switched-off limit with a NaN bit pattern (`no_limiter`); the twin uses `nothing`.
"""
struct HydrometeorProperties
    density::Float64
    mass_exponent::Float64
    mass_coefficient::Float64
    slope_bounds::Union{Nothing, NTuple{2, Float64}}
    minimum_mean_mass::Union{Nothing, Float64}
end

function HydrometeorProperties(
        density, mass_exponent; slope_bounds = nothing, minimum_mean_mass = nothing
)
    mass_coefficient = density * π * gamma_function(mass_exponent + 1.0) / 6.0
    return HydrometeorProperties(
        density, mass_exponent, mass_coefficient, slope_bounds, minimum_mean_mass
    )
end

"""
    hydrometeor_properties(ice_to_snow_size)

Properties of the six categories, as set by `micro_pumas_utils_init`.
`ice_to_snow_size` is the cloud ice to snow autoconversion diameter (m, `dcs`).
"""
function hydrometeor_properties(ice_to_snow_size)
    d = sphere_mass_exponent
    # 🐛🧊 Fortran writes 1._r8/2._r8*dcs, which is dcs/2, not 1/(2 dcs):
    #    see FINDINGS.md. The twin keeps the Fortran value.
    ice_slope_bounds = (1.0 / 2.0 * ice_to_snow_size, 1.0 / 1.0e-6)
    return (
        cloud_liquid = HydrometeorProperties(
            bulk_density.liquid_water, d; minimum_mean_mass = minimum_mean_mass
        ),
        cloud_ice = HydrometeorProperties(
            bulk_density.cloud_ice, d;
            slope_bounds = ice_slope_bounds, minimum_mean_mass = minimum_mean_mass
        ),
        rain = HydrometeorProperties(
            bulk_density.liquid_water, d; slope_bounds = rain_slope_bounds
        ),
        snow = HydrometeorProperties(
            bulk_density.snow, d; slope_bounds = snow_slope_bounds
        ),
        graupel = HydrometeorProperties(
            bulk_density.graupel, d; slope_bounds = snow_slope_bounds
        ),
        hail = HydrometeorProperties(
            bulk_density.hail, d; slope_bounds = snow_slope_bounds
        )
    )
end
