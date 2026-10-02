# 📏💧 Effective radius of rain and cloud liquid together.

"""
    add_rain_cloud_effective_radius(radius_sum, λ_rain, N₀_rain, λ_cloud, μ_cloud,
                                    q_rain, q_cloud, N_cloud)

Adds to `radius_sum` the effective radius (m) of rain plus cloud droplets,
r = 3 (q_r + q_c) / (4 ρ_w A), with the projected area A of rain
(N₀ π / (2 λ³)) plus droplets (N π (μ+1)(μ+2) / (4 λ²)); unchanged when A = 0.
PUMAS sums this over sub-steps and divides by their number later.
Fortran twin: `calc_rercld` in `micro_pumas_v1.F90` (`rercld`).
"""
function add_rain_cloud_effective_radius(radius_sum, λ_rain, N₀_rain, λ_cloud, μ_cloud,
        q_rain, q_cloud, N_cloud)
    area = λ_rain > 0.0 ? N₀_rain * π / (2.0 * λ_rain^3.0) : 0.0
    if λ_cloud > 0.0
        moment = rising_factorial(μ_cloud + 1.0, 2)
        area = area + N_cloud * π * moment / (4.0 * λ_cloud^2.0)
    end
    area > 0.0 || return radius_sum
    return radius_sum + 3.0 * (q_rain + q_cloud) / (4.0 * bulk_density.liquid_water * area)
end
