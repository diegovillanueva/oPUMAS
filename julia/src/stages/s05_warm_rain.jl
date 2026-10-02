# 🌧️🔁 Stage 5, autoconversion of cloud liquid to rain (micro_pumas_v1.F90 lines
#    2028–2165). KK2000 is always computed as a diagnostic.

"""
    autoconversion!(work, parameters)

🌧️ Stage 5. Autoconversion by `kk2000`, `sb2001` or the TAU bin code, chosen by
`warm_rain`.
"""
function autoconversion!(work, parameters)
    (; settings) = parameters
    (; inputs, state, incloud, rates, process_rates) = work
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    qr_precip = incloud.rain.q  # qric
    ρ = state.air.density  # rho
    relative_variance = inputs.relative_variance  # relvar
    μc = incloud.cloud_liquid.μ  # pgam
    qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc
    nr_autoconversion = rates.rain.autoconversion_from_cloud_liquid.number  # nprc
    nc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.number  # nprc1
    qc_kk2000 = process_rates.cloud_liquid.warm_rain_kk2000.mass  # qctend_KK2000
    nr_kk2000 = process_rates.rain.warm_rain_kk2000.number  # nrtend_KK2000
    nc_kk2000 = process_rates.cloud_liquid.warm_rain_kk2000.number  # nctend_KK2000
    qr_kk2000 = process_rates.rain.warm_rain_kk2000.mass  # qrtend_KK2000
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        kk = kk2000_autoconversion(parameters, qc_cloud[i, k], nc_cloud[i, k], ρ[i, k],
            relative_variance[i, k])
        qc_kk2000[i, k] = kk.mass
        nr_kk2000[i, k] = kk.rain_number
        nc_kk2000[i, k] = kk.droplet_number
    end
    if settings.warm_rain == :kk2000
        for k in 1:levels, i in 1:columns
            qc_autoconversion[i, k] = qc_kk2000[i, k]
            nc_autoconversion[i, k] = nc_kk2000[i, k]
            nr_autoconversion[i, k] = nr_kk2000[i, k]
            qr_kk2000[i, k] = -qc_kk2000[i, k]
        end
    elseif settings.warm_rain == :sb2001
        qc_sb2001 = process_rates.cloud_liquid.warm_rain_sb2001.mass  # qctend_SB2001
        nr_sb2001 = process_rates.rain.warm_rain_sb2001.number  # nrtend_SB2001
        nc_sb2001 = process_rates.cloud_liquid.warm_rain_sb2001.number  # nctend_SB2001
        qr_sb2001 = process_rates.rain.warm_rain_sb2001.mass  # qrtend_SB2001
        for k in 1:levels, i in 1:columns
            sb = sb2001_autoconversion(
                μc[i, k], qc_cloud[i, k], nc_cloud[i, k], qr_precip[i, k],
                ρ[i, k])
            qc_sb2001[i, k] = sb.mass
            nr_sb2001[i, k] = sb.rain_number
            nc_sb2001[i, k] = sb.droplet_number
        end
        for k in 1:levels, i in 1:columns
            qc_autoconversion[i, k] = qc_sb2001[i, k]
            nr_autoconversion[i, k] = nr_sb2001[i, k]
            nc_autoconversion[i, k] = nc_sb2001[i, k]
            qr_sb2001[i, k] = -qc_sb2001[i, k]
        end
    end
    settings.warm_rain == :tau && tau_warm_rain!(work, parameters)
    # 🧊🎛️ Sppe: scale the liquid autoconversion of kk2000 and sb2001
    factor = settings.autoconversion_factor
    if factor != 1.0 && settings.warm_rain in (:kk2000, :sb2001)
        for k in 1:levels, i in 1:columns
            qc_autoconversion[i, k] = qc_autoconversion[i, k] * factor
            nr_autoconversion[i, k] = nr_autoconversion[i, k] * factor
            nc_autoconversion[i, k] = nc_autoconversion[i, k] * factor
        end
    end
    return nothing
end

"""
    tau_warm_rain!(work, parameters)

🌧️🧮 All warm rain (autoconversion, accretion, self-collection) by the TAU bin code
([`tau_collection`](@ref)). Cloud losses become `prc` and `nprc1`; a rain number gain
becomes `nprc`, a loss (self-collection) `nragg`.
"""
function tau_warm_rain!(work, parameters)
    (; inputs, state, incloud, rates, process_rates) = work
    bins = parameters.tau_bins
    scratch = tau_scratch()
    Δt = inputs.time_step  # deltatin
    qc_grid = inputs.cloud_liquid.q  # qcn
    qr_grid = inputs.rain.q  # qrn
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    qr_precip = incloud.rain.q  # qric
    nr_precip = incloud.rain.N  # nric
    μc = incloud.cloud_liquid.μ  # pgam
    λc = incloud.cloud_liquid.λ  # lamc
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    liquid_fraction = incloud.liquid_fraction  # lcldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    ρ = state.air.density  # rho
    qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc
    nr_autoconversion = rates.rain.autoconversion_from_cloud_liquid.number  # nprc
    nc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.number  # nprc1
    nr_self_collection = rates.rain.self_collection.number  # nragg
    cloud = process_rates.cloud_liquid.warm_rain_tau
    rain = process_rates.rain.warm_rain_tau
    net_bin_gain = process_rates.liquid.warm_rain_tau.net_bin_gain  # gmnnn_lmnnn_TAU
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        tau = tau_collection(scratch, bins, Δt, ρ[i, k], qc_grid[i, k], qr_grid[i, k],
            qc_cloud[i, k], nc_cloud[i, k], qr_precip[i, k], nr_precip[i, k], μc[i, k],
            λc[i, k], λr[i, k], N₀r[i, k], liquid_fraction[i, k], precip_fraction[i, k])
        cloud.mass[i, k] = tau.tendency.cloud_mass  # qctend_TAU
        cloud.number[i, k] = tau.tendency.cloud_number  # nctend_TAU
        rain.mass[i, k] = tau.tendency.rain_mass  # qrtend_TAU
        rain.number[i, k] = tau.tendency.rain_number  # nrtend_TAU
        cloud.mass_after[i, k] = tau.after.cloud_mass  # qc_out_TAU
        cloud.number_after[i, k] = tau.after.cloud_number  # nc_out_TAU
        rain.mass_after[i, k] = tau.after.rain_mass  # qr_out_TAU
        rain.number_after[i, k] = tau.after.rain_number  # nr_out_TAU
        cloud.mass_scale[i, k] = tau.scale.cloud_mass  # scale_qc
        cloud.number_scale[i, k] = tau.scale.cloud_number  # scale_nc
        rain.mass_scale[i, k] = tau.scale.rain_mass  # scale_qr
        rain.number_scale[i, k] = tau.scale.rain_number  # scale_nr
        net_bin_gain[i, k] = tau.net_bin_gain
        cloud.mass_before[i, k] = qc_cloud[i, k]  # qc_in_TAU
        cloud.number_before[i, k] = nc_cloud[i, k]  # nc_in_TAU
        rain.mass_before[i, k] = qr_precip[i, k]  # qr_in_TAU
        rain.number_before[i, k] = nr_precip[i, k]  # nr_in_TAU
        qc_autoconversion[i, k] = -cloud.mass[i, k]
        nc_autoconversion[i, k] = -cloud.number[i, k]
        if rain.number[i, k] > 0.0
            nr_autoconversion[i, k] = rain.number[i, k]
        else
            nr_self_collection[i, k] = rain.number[i, k]
        end
    end
    return nothing
end
