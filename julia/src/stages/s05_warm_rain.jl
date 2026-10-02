# 🌧️🔁 Stage 5, autoconversion of cloud liquid to rain (micro_pumas_v1.F90 lines
#    2028–2165). KK2000 is always computed as a diagnostic.

"""
    autoconversion!(work, parameters)

🌧️ Stage 5. Autoconversion by `kk2000` or `sb2001`, chosen by `warm_rain`.
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
            sb = sb2001_autoconversion(μc[i, k], qc_cloud[i, k], nc_cloud[i, k], qr_precip[i, k],
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
    return nothing
end
