# 📏🌧️ Gamma size distributions N(D) = N₀ D^μ exp(-λ D) and helpers.

"""
    rising_factorial(x, n)

Rising factorial x (x+1) ⋯ (x+n-1) for integer `n`, and Γ(x+n)/Γ(x) for real `n`.
Fortran twin: `rising_factorial` in `micro_pumas_utils.F90`.
"""
function rising_factorial(x::Float64, n::Integer)
    product, factor = 1.0, x
    for _ in 1:n
        product = product * factor
        factor = factor + 1.0
    end
    return product
end
rising_factorial(x::Float64, n::Float64) = gamma_function(x + n) / gamma_function(x)

"""
    variance_coefficient(relative_variance, exponent)

Enhancement of a rate ∝ qᵃ by sub-grid variability of cloud water, for a gamma
distribution of q with inverse relative variance ν: Γ(ν+a) / (Γ(ν) νᵃ).
Fortran twin: `var_coef` in `micro_pumas_utils.F90`.
"""
function variance_coefficient(ν::Float64, a::Real)
    return rising_factorial(ν, a) / ν^a
end

"""
    latent_heat_factor(parameters, T, q_sat, L)

Psychrometric correction 1 + L² q_sat / (R_v T² c_p) of growth and evaporation rates
(dimensionless), for latent heat `L` (J kg⁻¹) and saturation specific humidity `q_sat`.
Fortran twin: `calc_ab` in `micro_pumas_utils.F90`.
"""
function latent_heat_factor(parameters, T, q_sat, L)
    (; gas_constant_water_vapor, specific_heat_dry_air) = parameters.host
    dqsat_dT = L * q_sat / (gas_constant_water_vapor * T^2)
    return 1.0 + dqsat_dT * L / specific_heat_dry_air
end

"""
    size_distribution(h, q, N) -> (; slope, number, intercept)

Exponential size distribution of a hydrometeor with properties `h`, in-cloud mass mixing
ratio `q` (kg kg⁻¹) and number `N` (kg⁻¹). Returns the slope λ (m⁻¹), the number after
the limits on mean mass and on λ, and the intercept N₀ = N λ. λ = 0 when q is negligible.

Fortran twin: `size_dist_param_basic` in `micro_pumas_utils.F90`.
"""
function size_distribution(h::HydrometeorProperties, q, N)
    if q <= smallest_mixing_ratio
        λ = 0.0
        return (; slope = λ, number = N, intercept = N * λ)
    end
    return limited_size_distribution(h, q, N, h.mass_coefficient, h.slope_bounds)
end

"""
Slope λ = (c N / q)^(1/d) clamped to `slope_bounds`, with N adjusted to the clamp.
Shared by [`size_distribution`](@ref) and [`droplet_size_distribution`](@ref).
"""
function limited_size_distribution(h::HydrometeorProperties, q, N, c, slope_bounds)
    d = h.mass_exponent
    if h.minimum_mean_mass !== nothing
        N = min(N, q / h.minimum_mean_mass)
    end
    λ = (c * N / q)^(1.0 / d)
    if slope_bounds !== nothing
        λ_min, λ_max = slope_bounds
        if λ < λ_min
            λ = λ_min
            N = λ^d * q / c
        elseif λ > λ_max
            λ = λ_max
            N = λ^d * q / c
        end
    end
    return (; slope = λ, number = N, intercept = N * λ)
end

"""
    droplet_size_distribution(h, q, N, ρ) -> (; shape, slope, number)

Gamma size distribution of cloud droplets: shape μ from the fit of Rotstayn and Liu
(2003) to droplet number, then slope λ with mean diameter limited to 2–50 μm.
`ρ` is the air density (kg m⁻³). When q is negligible, μ = -100 (unusable) and λ = 0.

Fortran twin: `size_dist_param_liq` in `micro_pumas_utils.F90`.
"""
function droplet_size_distribution(h::HydrometeorProperties, q, N, ρ)
    q > smallest_mixing_ratio || return (; shape = -100.0, slope = 0.0, number = N)
    dispersion = 1.0 - 0.7 * exp(-0.008 * 1.0e-6 * N * ρ)
    μ = max(1.0 / (dispersion^2) - 1.0, 2.0)
    c = π / 6.0 * h.density * droplet_moment_factor(h, μ + 1.0)
    bounds = ((μ + 1.0) * 1.0 / 50.0e-6, (μ + 1.0) * 1.0 / 2.0e-6)
    sd = limited_size_distribution(h, q, N, c, bounds)
    return (; shape = μ, slope = sd.slope, number = sd.number)
end

"Γ(μ+1+d)/Γ(μ+1), the integer product when the mass exponent d is 3."
function droplet_moment_factor(h::HydrometeorProperties, μ_plus_1)
    h.mass_exponent == 3.0 && return rising_factorial(μ_plus_1, 3)
    return rising_factorial(μ_plus_1, h.mass_exponent)
end

"""
    mean_diameter(q, N, ρ_air, ρ_particle)

Mean diameter (m) of an exponential distribution of spheres of density `ρ_particle`, from
mass `q` (kg kg⁻¹) and number `N`. Fortran twin: `avg_diameter` in `micro_pumas_utils.F90`.
"""
mean_diameter(q, N, ρ_air, ρ_particle) = (q * ρ_air / (π * ρ_particle * N))^(1.0 / 3.0)
