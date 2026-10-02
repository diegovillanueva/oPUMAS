# 🌧️🧮 TAU bin warm rain (pumas_stochastic_collect_tau.F90, Gettelman and Chen 2018):
#    cloud and rain size distributions are cut into 35 mass-doubling bins, collide for one
#    time step with a stochastic collection kernel, and are summed back into cloud and
#    rain. Bins use the Fortran's cgs units (g, cm) unless the name says otherwise.

"Number of drop size bins (`ncd`)."
const tau_bin_count = 35

"""
Mass of a 6 μm drop, kg, given to drops that have mass but no number (`m1`). Written as
the value the Fortran compiler folds (FINDINGS.md), one unit in the last place below the
run-time value of 4/3 π ρ_w (6 μm)³.
"""
const tau_new_drop_mass = 9.047786842338604e-13

"Substeps of the collision loop within one time step (`sub_step`)."
const tau_substeps = 60

"""
    tau_bins(kernel_file) -> NamedTuple

Mass-doubling bin grid of Graham Feingold and the collection kernel read from
`kernel_file`:

- `edge_mass` (g) and `edge_diameter` (cm): the 36 bin edges, from a 3.125 μm drop
- `mean_mass` (g) and `mean_diameter` (cm): the 35 bin mid-points
- `cutoff`: first edge at or above 40 μm, the cloud-rain boundary (`cutoff_id`)
- `kernel`: K(i, j) (lower triangle) times 1.5 (m_i + m_j), negative values set to 0

Fortran twin: `pumas_stochastic_kernel_init`, which calls `calc_bins` (`medge`,
`diamedge`, `mmean`, `diammean`, `cutoff_id`, `knn`).
"""
function tau_bins(kernel_file)
    ρ_w = bulk_density.liquid_water
    n = tau_bin_count
    edge_diameter = zeros(n + 1)  # DIAM
    edge_mass = zeros(n + 1)  # X, radsl
    edge_diameter[1] = 1.5625 * 2.0e-4
    edge_mass[1] = π / 6.0 * edge_diameter[1]^3 * ρ_w / 1000.0
    for l in 2:(n + 1)
        edge_mass[l] = 2.0 * edge_mass[l - 1]
        edge_diameter[l] = (6.0 / π * edge_mass[l] * 1000.0 / ρ_w)^(1.0 / 3.0)
    end
    mean_mass = zeros(n)  # radl, mmean
    mean_diameter = zeros(n)  # diammean
    for l in 1:n
        mean_mass[l] = (edge_mass[l] + edge_mass[l + 1]) / 2.0
        mean_diameter[l] = (6.0 / π * mean_mass[l] * 1000.0 / ρ_w)^(1.0 / 3.0)
    end
    cutoff = 0
    for l in (n + 1):-1:1
        edge_diameter[l] >= 40.0e-4 && (cutoff = l)
    end
    # 📥 one value per line in columns 3–14 (FORMAT(2X,E12.5)), lower triangle by rows
    kernel = zeros(n, n)  # knn
    lines = Iterators.Stateful(eachline(kernel_file))
    for i in 1:n, j in 1:i
        kernel[i, j] = parse(Float64, strip(popfirst!(lines)[3:14]))
        kernel[i, j] = (mean_mass[i] * 1.5 + mean_mass[j] * 1.5) * kernel[i, j]
        kernel[i, j] < 0.0 && (kernel[i, j] = 0.0)
    end
    return (; edge_mass, edge_diameter, mean_mass, mean_diameter, cutoff, kernel)
end

"Work vectors of one grid point for [`tau_collection`](@ref), allocated once."
function tau_scratch()
    bins() = zeros(tau_bin_count)
    return (
        cloud = (mass = bins(), number = bins()),  # amk_c ank_c
        rain = (mass = bins(), number = bins()),  # amk_r ank_r
        before = (mass = bins(), number = bins()),  # amk ank
        after = (mass = bins(), number = bins()),  # amk0 ank0
        gain = (mass = bins(), number = bins()),  # gmnnn gnnnn
        loss = (mass = bins(), number = bins()),  # lmnnn lnnnn
        step_gain = (mass = bins(), number = bins()),  # gmnnn0 gnnnn0
        step_loss = (mass = bins(), number = bins()),  # lmnnn0 lnnnn0
        # 📐 smn0 snn0 smcn, amn2 amn3, psin fn, fpsin xpsin hpsin fn2
        moments = (mass = bins(), number = bins(), mean_mass = bins(),
            second = bins(), third = bins(), slope = bins(), offset = bins(),
            fpsin = bins(), xpsin = bins(), hpsin = bins(), half_offset = bins())
    )
end

"""
    tau_bin_distribution!(scratch, bins, q_cloud_grid, q_rain_grid, q_cloud, N_cloud,
        q_rain, N_rain, μ_cloud, λ_cloud, λ_rain, N₀_rain, liquid_fraction,
        precip_fraction)

Fill `scratch.cloud`, `scratch.rain` and `scratch.before` (kg kg⁻¹ and kg⁻¹ per bin)
from the gamma droplet and exponential rain distributions, each rescaled so that the
bins hold exactly the in-cloud mass and number, times the cloud or precipitation
fraction. Bins are only filled where the grid-mean and in-cloud mass exceed
`smallest_mixing_ratio`.

Returns `(; scale, cutoff)`: the four rescaling factors (bin sum / moment) and the bin
of least mass between the cloud and rain peaks (0 when there is no such minimum).
Fortran twin: `cam_bin_distribute` (`scale_qc`, `scale_nc`, `scale_qr`, `scale_nr`,
`cutoff_amk`).
"""
function tau_bin_distribution!(scratch, bins, q_cloud_grid, q_rain_grid, q_cloud,
        N_cloud, q_rain, N_rain, μ_cloud, λ_cloud, λ_rain, N₀_rain, liquid_fraction,
        precip_fraction)
    (; cloud, rain, before) = scratch
    (; mean_diameter, edge_diameter, mean_mass) = bins
    for j in 1:tau_bin_count
        cloud.number[j] = 0.0
        cloud.mass[j] = 0.0
        rain.number[j] = 0.0
        rain.mass[j] = 0.0
    end
    cloud_number_scale = 0.0  # scale_nc
    cloud_mass_scale = 0.0  # scale_qc
    rain_number_scale = 0.0  # scale_nr
    rain_mass_scale = 0.0  # scale_qr
    cutoff = 0
    cloud_peak = 0  # id_max_qc
    rain_peak = 0  # id_max_qr
    cloud_peak_mass = 0.0  # max_qc
    rain_peak_mass = 0.0  # max_qr
    if q_cloud_grid > smallest_mixing_ratio && q_cloud > smallest_mixing_ratio
        for j in 1:tau_bin_count
            density = N_cloud * λ_cloud^(μ_cloud + 1.0) / gamma_function(μ_cloud + 1.0) *
                      (mean_diameter[j] * 1.0e-2)^μ_cloud *
                      exp(-λ_cloud * mean_diameter[j] * 1.0e-2)  # phi, D cm → m
            cloud.number[j] = density * (edge_diameter[j + 1] - edge_diameter[j]) * 1.0e-2
            cloud.mass[j] = density * (edge_diameter[j + 1] - edge_diameter[j]) * 1.0e-2 *
                            mean_mass[j] * 1.0e-3
            cloud_number_scale = cloud_number_scale + cloud.number[j]
            cloud_mass_scale = cloud_mass_scale + cloud.mass[j]
        end
        cloud_number_scale = cloud_number_scale / N_cloud
        cloud_mass_scale = cloud_mass_scale / q_cloud
        for j in 1:tau_bin_count
            cloud.number[j] = cloud.number[j] / cloud_number_scale * liquid_fraction
            cloud.mass[j] = cloud.mass[j] / cloud_mass_scale * liquid_fraction
            if cloud.mass[j] > cloud_peak_mass
                cloud_peak = j
                cloud_peak_mass = cloud.mass[j]
            end
        end
    end
    if q_rain_grid > smallest_mixing_ratio && q_rain > smallest_mixing_ratio
        for j in 1:tau_bin_count
            density = N₀_rain * exp(-λ_rain * mean_diameter[j] * 1.0e-2)  # phi
            rain.number[j] = density * (edge_diameter[j + 1] - edge_diameter[j]) * 1.0e-2
            rain.mass[j] = density * (edge_diameter[j + 1] - edge_diameter[j]) * 1.0e-2 *
                           mean_mass[j] * 1.0e-3
            rain_number_scale = rain_number_scale + rain.number[j]
            rain_mass_scale = rain_mass_scale + rain.mass[j]
        end
        rain_number_scale = rain_number_scale / N_rain
        rain_mass_scale = rain_mass_scale / q_rain
        for j in 1:tau_bin_count
            rain.number[j] = rain.number[j] / rain_number_scale * precip_fraction
            rain.mass[j] = rain.mass[j] / rain_mass_scale * precip_fraction
            if rain.mass[j] > rain_peak_mass
                rain_peak = j
                rain_peak_mass = rain.mass[j]
            end
        end
    end
    for j in 1:tau_bin_count
        before.mass[j] = cloud.mass[j] + rain.mass[j]
        before.number[j] = cloud.number[j] + rain.number[j]
    end
    if cloud_peak > 0 && rain_peak > 0
        ratio = cloud_peak_mass / rain_peak_mass
        if ratio < 10.0 || ratio > 0.1  # as in the Fortran: always true
            least = before.mass[cloud_peak]  # min_amk
            for j in cloud_peak:rain_peak
                if before.mass[j] <= least
                    cutoff = j
                    least = before.mass[j]
                end
            end
        end
    end
    scale = (cloud_mass = cloud_mass_scale, cloud_number = cloud_number_scale,
        rain_mass = rain_mass_scale, rain_number = rain_number_scale)
    return (; scale, cutoff)
end

"""
    tau_collision_rates!(scratch, bins, ρ_air)

Gain and loss of mass (g cm⁻³ s⁻¹) and number (cm⁻³ s⁻¹) of each bin by drop-drop
collection, for the bin contents `scratch.after` (kg kg⁻¹, kg⁻¹), written to
`scratch.step_gain` and `scratch.step_loss`. Within a bin, number is linear in mass
between the edges (method of moments of Tzivion et al. 1987).
Fortran twin: `compute_coll_params` (`gmnnn`, `gnnnn`, `lmnnn`, `lnnnn`).
"""
function tau_collision_rates!(scratch, bins, ρ_air)
    m = scratch.moments
    q, N = scratch.after.mass, scratch.after.number
    gain, loss = scratch.step_gain, scratch.step_loss
    x = bins.edge_mass  # xk_gr
    K = bins.kernel  # knn
    tiny = 1.0e-12  # lsmall
    n = tau_bin_count
    for l in 1:n
        m.mass[l] = q[l] * ρ_air / 1.0e3  # smn0
        m.number[l] = N[l] * ρ_air / 1.0e6  # snn0
        if m.mass[l] < tiny || m.number[l] < tiny
            m.mass[l] = 0.0
            m.number[l] = 0.0
        end
    end
    for l in 1:n
        if m.mass[l] > 0.0 && m.number[l] > 0.0
            m.mean_mass[l] = m.mass[l] / m.number[l]
            m.mean_mass[l] > 2.0 * x[l] && (m.mean_mass[l] = 2.0 * x[l])
            m.mean_mass[l] < x[l] && (m.mean_mass[l] = x[l])
        else
            m.mean_mass[l] = 0.0
        end
        if m.mean_mass[l] < x[l] || m.mean_mass[l] > 2.0 * x[l] || m.mean_mass[l] == 0.0
            shape = 1.0  # apn
        else
            shape = 0.5 * (1.0 + 3.0 * (x[l] / m.mean_mass[l]) -
                     2.0 * ((x[l] / m.mean_mass[l])^2.0))
        end
        if m.number[l] > tiny
            m.second[l] = shape * m.mass[l] * m.mass[l] / m.number[l]
            m.third[l] = shape * shape * shape * m.mass[l] * m.mass[l] * m.mass[l] /
                         (m.number[l] * m.number[l])
        else
            m.second[l] = 0.0
            m.third[l] = 0.0
        end
        if m.mean_mass[l] < x[l]
            m.slope[l] = 0.0
            m.offset[l] = 2.0 * m.number[l] / x[l]
        elseif m.mean_mass[l] > 2.0 * x[l]
            m.offset[l] = 0.0
            m.slope[l] = 2.0 * m.number[l] / x[l]
        else
            m.slope[l] = 2.0 / x[l] * (m.mass[l] / x[l] - m.number[l])
            m.offset[l] = 2.0 / x[l] * (2.0 * m.number[l] - m.mass[l] / x[l])
        end
        if m.number[l] < tiny || m.mass[l] < tiny
            m.slope[l] = 0.0
            m.offset[l] = 0.0
        end
        m.fpsin[l] = 0.5 / x[l] * (m.slope[l] - m.offset[l])
        m.xpsin[l] = 2.0 * x[l] * m.slope[l]
        m.hpsin[l] = m.slope[l] - 0.5 * m.offset[l]
        m.half_offset[l] = m.offset[l] / 2.0
    end
    psin, fpsin, xpsin, hpsin, fn2 = m.slope, m.fpsin, m.xpsin, m.hpsin, m.half_offset
    smn0, snn0, amn2, amn3 = m.mass, m.number, m.second, m.third
    fill!(gain.mass, 0.0)
    fill!(gain.number, 0.0)
    fill!(loss.mass, 0.0)
    fill!(loss.number, 0.0)
    # 🔁 drops of bin l-1 collecting smaller drops grow into bin l
    for l in 3:(n - 1)
        lm1 = l - 1
        for ll in 1:(l - 2)
            gain.number[l] = gain.number[l] +
                             (psin[lm1] * smn0[ll] - fpsin[lm1] * amn2[ll]) * K[lm1, ll]
            gain.mass[l] = gain.mass[l] +
                           (x[l] * psin[lm1] * smn0[ll] + fn2[lm1] * amn2[ll] -
                            fpsin[lm1] * amn3[ll]) * K[lm1, ll]
        end
    end
    for l in 2:(n - 1)
        lm1 = l - 1
        gain.number[l] = gain.number[l] + 0.5 * snn0[lm1] * snn0[lm1] * K[lm1, lm1]
        gain.mass[l] = gain.mass[l] +
                       0.5 * (snn0[lm1] * smn0[lm1] + smn0[lm1] * snn0[lm1]) * K[lm1, lm1]
        for ll in 1:(l - 1)
            loss.number[l] = loss.number[l] +
                             (psin[l] * smn0[ll] - fpsin[l] * amn2[ll]) * K[l, ll]
            gain.mass[l] = gain.mass[l] +
                           (smn0[ll] * snn0[l] - psin[l] * amn2[ll] +
                            fpsin[l] * amn3[ll]) * K[l, ll]
            loss.mass[l] = loss.mass[l] +
                           (xpsin[l] * smn0[ll] - hpsin[l] * amn2[ll]) * K[l, ll]
        end
    end
    for l in 1:(n - 1), ll in l:(n - 1)
        loss.number[l] = loss.number[l] + (snn0[ll] * snn0[l]) * K[ll, l]
        loss.mass[l] = loss.mass[l] + (snn0[ll] * smn0[l]) * K[ll, l]
    end
    return nothing
end

"""
    tau_collection(scratch, bins, Δt, ρ_air, q_cloud_grid, q_rain_grid, q_cloud,
        N_cloud, q_rain, N_rain, μ_cloud, λ_cloud, λ_rain, N₀_rain, liquid_fraction,
        precip_fraction)

🌧️ Autoconversion, accretion and self-collection of liquid drops by the TAU bin code for
one grid point and one time step Δt (s): bin the distributions
([`tau_bin_distribution!`](@ref)), let them collide in 60 substeps
([`tau_collision_rates!`](@ref)) with the loss rescaled to the gain so that each substep
conserves mass, and sum the bins below and above the cut back into cloud and rain. The
cut is the bin of least mass between the cloud and rain peaks, else 40 μm. The result is
clipped so that cloud only loses, rain only gains mass, and cloud plus rain mass is
conserved; a category with mass but no number gets 6 μm drops.

Returns in-cloud values: `tendency` (cloud and rain mass, kg kg⁻¹ s⁻¹, and number,
kg⁻¹ s⁻¹; cloud ≤ 0, rain mass ≥ 0), `after` (the four values after the step), `scale`
(see [`tau_bin_distribution!`](@ref)) and `net_bin_gain` (sum over bins of gain - loss,
g cm⁻³ s⁻¹).
Fortran twin: `pumas_stochastic_collect_tau_tend` (`qctend_TAU`, `nctend_TAU`,
`qrtend_TAU`, `nrtend_TAU`, `qcin_new`, `ncin_new`, `qrin_new`, `nrin_new`,
`gmnnn_lmnnn_TAU`).
"""
function tau_collection(scratch, bins, Δt, ρ_air, q_cloud_grid, q_rain_grid, q_cloud,
        N_cloud, q_rain, N_rain, μ_cloud, λ_cloud, λ_rain, N₀_rain, liquid_fraction,
        precip_fraction)
    (; before, after, gain, loss, step_gain, step_loss) = scratch
    n = tau_bin_count
    cut = bins.cutoff - 1
    distribution = tau_bin_distribution!(scratch, bins, q_cloud_grid, q_rain_grid,
        q_cloud, N_cloud, q_rain, N_rain, μ_cloud, λ_cloud, λ_rain, N₀_rain,
        liquid_fraction, precip_fraction)
    distribution.cutoff > 0 && (cut = distribution.cutoff)
    for field in (gain.mass, gain.number, loss.mass, loss.number)
        fill!(field, 0.0)
    end
    copyto!(after.mass, before.mass)
    copyto!(after.number, before.number)
    for _ in 1:tau_substeps
        tau_collision_rates!(scratch, bins, ρ_air)
        total_gain = 0.0  # all_gmnnn
        total_loss = 0.0  # all_lmnnn
        for l in 1:n
            total_gain = total_gain + step_gain.mass[l]
            total_loss = total_loss + step_loss.mass[l]
        end
        if total_gain == 0.0 || total_loss == 0.0
            fill!(step_gain.mass, 0.0)
            fill!(step_loss.mass, 0.0)
        else
            for l in 1:n
                step_loss.mass[l] = step_loss.mass[l] * (total_gain / total_loss)
            end
        end
        for l in 1:n
            after.mass[l] = after.mass[l] +
                            (step_gain.mass[l] - step_loss.mass[l]) * 1.0e3 / ρ_air * Δt /
                            Float64(tau_substeps)
            after.number[l] = after.number[l] +
                              (step_gain.number[l] - step_loss.number[l]) * 1.0e6 / ρ_air *
                              Δt / Float64(tau_substeps)
            gain.mass[l] = gain.mass[l] + step_gain.mass[l] / tau_substeps
            gain.number[l] = gain.number[l] + step_gain.number[l] / tau_substeps
            loss.mass[l] = loss.mass[l] + step_loss.mass[l] / tau_substeps
            loss.number[l] = loss.number[l] + step_loss.number[l] / tau_substeps
        end
    end
    qc_old = nc_old = qr_old = nr_old = 0.0
    qc_new = nc_new = qr_new = nr_new = 0.0
    net_bin_gain = 0.0  # gmnnn_lmnnn_TAU
    for l in 1:cut  # 💧 cloud bins
        qc_old = qc_old + before.mass[l]
        nc_old = nc_old + before.number[l]
        qc_new = qc_new + (gain.mass[l] - loss.mass[l]) * 1.0e3 / ρ_air * Δt
        nc_new = nc_new + (gain.number[l] - loss.number[l]) * 1.0e6 / ρ_air * Δt
        net_bin_gain = net_bin_gain + gain.mass[l] - loss.mass[l]
    end
    for l in (cut + 1):n  # 🌧️ rain bins
        qr_old = qr_old + before.mass[l]
        nr_old = nr_old + before.number[l]
        qr_new = qr_new + (gain.mass[l] - loss.mass[l]) * 1.0e3 / ρ_air * Δt
        nr_new = nr_new + (gain.number[l] - loss.number[l]) * 1.0e6 / ρ_air * Δt
        net_bin_gain = net_bin_gain + gain.mass[l] - loss.mass[l]
    end
    qc_new = qc_new + qc_old
    nc_new = nc_new + nc_old
    qr_new = qr_new + qr_old
    nr_new = nr_new + nr_old
    # ⚖️ conservation: no negatives, cloud only loses, rain only gains, qc + qr kept
    qc_new = max(qc_new, 0.0)
    nc_new = max(nc_new, 0.0)
    qr_new = max(qr_new, 0.0)
    nr_new = max(nr_new, 0.0)
    qc_new = min(qc_new, q_cloud)
    nc_new = min(nc_new, N_cloud)
    qr_new = max(qr_new, q_rain)
    if qc_new + qr_new > 0.0
        mass_scale = (q_cloud + q_rain) / (qc_new + qr_new)  # qscl
    else
        mass_scale = 0.0
    end
    qc_new = qc_new * mass_scale
    qr_new = qr_new * mass_scale
    qc_new < smallest_mixing_ratio && (nc_new = 0.0)
    qr_new < smallest_mixing_ratio && (nr_new = 0.0)
    if qc_new > smallest_mixing_ratio && nc_new < smallest_mixing_ratio
        nc_new = qc_new / tau_new_drop_mass
    end
    if qr_new > smallest_mixing_ratio && nr_new < smallest_mixing_ratio
        nr_new = qr_new / tau_new_drop_mass
    end
    tendency = (
        cloud_mass = min((qc_new - q_cloud) / Δt, 0.0),
        cloud_number = min((nc_new - N_cloud) / Δt, 0.0),
        rain_mass = max((qr_new - q_rain) / Δt, 0.0),
        rain_number = (nr_new - N_rain) / Δt
    )
    after_step = (cloud_mass = qc_new, cloud_number = nc_new, rain_mass = qr_new,
        rain_number = nr_new)
    return (; tendency, after = after_step, distribution.scale, net_bin_gain)
end
