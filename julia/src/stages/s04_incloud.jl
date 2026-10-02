# ☁️🔍 Stage 4, in-cloud values and precipitation fraction (micro_pumas_v1.F90 lines
#    1899–2027), with the size distributions of cloud liquid and rain.

"""
    incloud_values!(work, parameters)

☁️ Stage 4. Grid-mean cloud liquid and ice become in-cloud values, rain becomes an
in-precipitation value, and the precipitation fraction follows the method of the
settings.
"""
function incloud_values!(work, parameters)
    (; settings, hydrometeors) = parameters
    (; state, incloud, process_rates) = work
    μc_out = process_rates.cloud_liquid.size_distribution.shape  # pgam_out
    N₀r_out = process_rates.rain.size_distribution.intercept  # n0r_out
    λc_out = process_rates.cloud_liquid.size_distribution.slope  # lamc_out
    λr_out = process_rates.rain.size_distribution.slope  # lamr_out
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    qi = state.cloud_ice.q
    ni = state.cloud_ice.N
    qr = state.rain.q
    nr = state.rain.N
    qs = state.snow.q
    qg = state.graupel.q
    ρ = state.air.density  # rho
    cloud_fraction = incloud.cloud_fraction  # cldm
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    qi_cloud = incloud.cloud_ice.q  # qiic
    ni_cloud = incloud.cloud_ice.N  # niic
    qr_precip = incloud.rain.q  # qric
    nr_precip = incloud.rain.N  # nric
    μc = incloud.cloud_liquid.μ  # pgam
    λc = incloud.cloud_liquid.λ  # lamc
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= qsmall
            qc_cloud[i, k] = min(qc[i, k] / liquid_fraction[i, k], 5.0e-3)
            nc_cloud[i, k] = max(nc[i, k] / liquid_fraction[i, k], 0.0)
            if settings.constant_droplet_number
                nc_cloud[i, k] = settings.droplet_number / ρ[i, k]
            end
        else
            qc_cloud[i, k] = 0.0
            nc_cloud[i, k] = 0.0
        end
        if qi[i, k] >= qsmall
            qi_cloud[i, k] = min(qi[i, k] / ice_fraction[i, k], 5.0e-3)
            ni_cloud[i, k] = max(ni[i, k] / ice_fraction[i, k], 0.0)
            if settings.constant_ice_number
                ni_cloud[i, k] = settings.ice_number / ρ[i, k]
            end
        else
            qi_cloud[i, k] = 0.0
            ni_cloud[i, k] = 0.0
        end
        precip_fraction[i, k] = cloud_fraction[i, k]
    end

    # 🌧️ precipitation fraction, top down
    if settings.precipitation_fraction_method == :in_cloud
        for i in 1:columns, k in 2:levels
            if qc[i, k] < qsmall && qi[i, k] < qsmall
                precip_fraction[i, k] = precip_fraction[i, k-1]
            end
        end
    else
        for i in 1:columns, k in 2:levels
            if qr[i, k-1] >= qsmall || qs[i, k-1] >= qsmall || qg[i, k-1] >= qsmall
                precip_fraction[i, k] = max(precip_fraction[i, k-1], precip_fraction[i, k])
            end
        end
    end

    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, qc_cloud[i, k],
            nc_cloud[i, k], ρ[i, k])
        μc[i, k], λc[i, k], nc_cloud[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        qr_precip[i, k] = qr[i, k] / precip_fraction[i, k]
        nr_precip[i, k] = nr[i, k] / precip_fraction[i, k]
        qr_precip[i, k] = min(qr_precip[i, k], 0.01)
        if qr_precip[i, k] < qsmall
            qr_precip[i, k] = 0.0
            nr_precip[i, k] = 0.0
        end
        nr_precip[i, k] = max(nr_precip[i, k], 0.0)
    end
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qr_precip[i, k], nr_precip[i, k])
        λr[i, k], N₀r[i, k], nr_precip[i, k] = rain.slope, rain.intercept, rain.number
    end
    for k in 1:levels, i in 1:columns
        μc_out[i, k] = μc[i, k]
        N₀r_out[i, k] = N₀r[i, k]
        λc_out[i, k] = λc[i, k]
        λr_out[i, k] = λr[i, k]
    end
    return nothing
end
