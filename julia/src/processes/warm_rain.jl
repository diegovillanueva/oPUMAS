# 🌧️💧 Warm rain: cloud droplets turning into rain, and rain drops colliding.
#    Rates are per second, for in-cloud (or in-precipitation) mixing ratios.

"""
    kk2000_autoconversion(parameters, q_cloud, N_cloud, ρ_air, relative_variance)

Autoconversion of cloud droplets to rain after Khairoutdinov and Kogan (2000):

    rate = E · f · 1350 · q_cᵃ · (N_c ρ_air 10⁻⁶)ᵇ

with E the sub-grid variance enhancement, f the prefactor, a and b the water and
number exponents of the settings. New raindrops have a 25 μm radius.

Returns `(; mass, rain_number, droplet_number)`: the cloud water converted
(kg kg⁻¹ s⁻¹), raindrops made and droplets lost (kg⁻¹ s⁻¹), all ≥ 0.
Nothing happens below `small_incloud_mixing_ratio`.
Fortran twin: `kk2000_liq_autoconversion` (`prc`, `nprc`, `nprc1`).
"""
function kk2000_autoconversion(parameters, q_cloud, N_cloud, ρ_air, relative_variance)
    s = parameters.settings
    q_cloud >= small_incloud_mixing_ratio ||
        return (; mass = 0.0, rain_number = 0.0, droplet_number = 0.0)
    a = s.autoconversion_water_exponent
    enhancement =
        s.uniform_subcolumns ? 1.0 : variance_coefficient(relative_variance, a)
    mass = enhancement * s.autoconversion_prefactor * 1350.0 * q_cloud^a *
           (N_cloud * 1.0e-6 * ρ_air)^s.autoconversion_number_exponent
    return (;
        mass,
        rain_number = mass * (1.0 / autoconversion_drop_mass_25um),
        droplet_number = mass * N_cloud / q_cloud,
    )
end

"""
Mass spectral shape ν of cloud droplets for each droplet size shape parameter μ =
1, 2, …, 16, from Seifert and Beheng (2001) (`dnu`). Linear in between.
"""
const droplet_mass_shape_table = (
    0.0, -0.557, -0.430, -0.307, -0.186, -0.067, 0.050, 0.167,
    0.282, 0.397, 0.512, 0.626, 0.739, 0.853, 0.966, 0.966,
)

"Autoconversion and accretion constants of Seifert and Beheng (2001) (`kc`, `kr`)."
const seifert_beheng = (autoconversion = 9.44e9, accretion = 5.78e3)

"""
    sb2001_autoconversion(shape, q_cloud, N_cloud, q_rain, ρ_air)

Autoconversion of cloud droplets to rain after Seifert and Beheng (2001), for droplet
size shape parameter `shape` (μ). The rate grows with the rain fraction of liquid water
through the universal function Φ(τ) = 600 τ^0.68 (1 - τ^0.68)³, τ = q_r / (q_c + q_r).
New raindrops have a 40 μm radius.

Returns `(; mass, rain_number, droplet_number)` as [`kk2000_autoconversion`](@ref).
Nothing happens at or below `smallest_mixing_ratio`.
Fortran twin: `sb2001v2_liq_autoconversion` (`au`, `nprc`, `nprc1`).
"""
function sb2001_autoconversion(shape, q_cloud, N_cloud, q_rain, ρ_air)
    q_cloud > smallest_mixing_ratio ||
        return (; mass = 0.0, rain_number = 0.0, droplet_number = 0.0)
    k = trunc(Int, shape)
    ν_table = droplet_mass_shape_table
    ν = ν_table[k] + (ν_table[k+1] - ν_table[k]) * (shape - k)
    # 🧮 τ: rain fraction of liquid water, written as the Fortran does.
    τ = 1.0 - q_cloud / (q_cloud + q_rain)
    Φ = 600.0 * τ^0.68 * (1.0 - τ^0.68)^3
    x_star = 2.6e-7   # ⚖️ drop mass x* separating cloud from rain, in g (cgs units)
    mass = seifert_beheng.autoconversion / (20.0 * x_star) *
           (ν + 2.0) * (ν + 4.0) / (ν + 1.0)^2.0 *
           (ρ_air * q_cloud / 1000.0)^4.0 / (ρ_air * N_cloud / 1.0e6)^2.0 *
           (1.0 + Φ / (1.0 - τ)^2) * 1000.0 / ρ_air
    return (;
        mass,
        rain_number = mass / autoconversion_drop_mass_40um,
        droplet_number = mass * 2.0 / x_star * 1000.0,
    )
end

"""
    sb2001_accretion(q_cloud, N_cloud, q_rain, ρ_air)

Accretion of cloud droplets by rain after Seifert and Beheng (2001), damped by
(τ / (τ + 5 10⁻⁴))⁴ while rain is a tiny fraction τ of liquid water.

Returns `(; mass, droplet_number)`: cloud water collected (kg kg⁻¹ s⁻¹) and droplets
lost (kg⁻¹ s⁻¹). Fortran twin: `sb2001v2_accre_cld_water_rain` (`pra`, `npra`).
"""
function sb2001_accretion(q_cloud, N_cloud, q_rain, ρ_air)
    q_cloud > smallest_mixing_ratio || return (; mass = 0.0, droplet_number = 0.0)
    τ = 1.0 - q_cloud / (q_cloud + q_rain)
    damping = (τ / (τ + 5.0e-4))^4.0
    mass = seifert_beheng.accretion * ρ_air * 0.001 * q_cloud * q_rain * damping
    droplet_number = mass * ρ_air * 0.001 * (N_cloud * ρ_air * 1.0e-6) /
                     (q_cloud * ρ_air * 0.001) * 1.0e6 / ρ_air
    return (; mass, droplet_number)
end

"""
    accretion(parameters, q_rain, q_cloud, N_cloud, relative_variance, enhancement)

Accretion of cloud droplets by rain after Khairoutdinov and Kogan (2000):

    rate = e · E · 67 (q_c q_r)^1.15

with `enhancement` e from the host and E the sub-grid variance enhancement (both 1 for
uniform sub-columns). Returns `(; mass, droplet_number)`: cloud water collected
(kg kg⁻¹ s⁻¹) and droplets lost (kg⁻¹ s⁻¹).
Fortran twin: `accrete_cloud_water_rain` (`pra`, `npra`).
"""
function accretion(parameters, q_rain, q_cloud, N_cloud, relative_variance, enhancement)
    if q_rain < smallest_mixing_ratio || q_cloud < smallest_mixing_ratio
        return (; mass = 0.0, droplet_number = 0.0)
    end
    factor = parameters.settings.uniform_subcolumns ? 1.0 :
             enhancement * variance_coefficient(relative_variance, 1.15)
    mass = factor * 67.0 * (q_cloud * q_rain)^1.15
    return (; mass, droplet_number = mass * N_cloud / q_cloud)
end

"""
    rain_self_collection(ρ_air, q_rain, N_rain)

Loss of raindrop number by raindrops merging, -8 N_r q_r ρ_air (kg⁻¹ s⁻¹), as in
Beheng (1994). Fortran twin: `self_collection_rain` (`nragg`).
"""
function rain_self_collection(ρ_air, q_rain, N_rain)
    q_rain >= smallest_mixing_ratio || return 0.0
    return -8.0 * N_rain * q_rain * ρ_air
end
