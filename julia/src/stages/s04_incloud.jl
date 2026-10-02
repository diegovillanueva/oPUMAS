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
    F = work.fortran
    (; qc, nc, qi, ni, qr, nr, qs, qg, rho, cldm, lcldm, icldm, precip_frac) = F
    (; qcic, ncic, qiic, niic, qric, nric, pgam, lamc, lamr, n0r) = F
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= qsmall
            qcic[i, k] = min(qc[i, k] / lcldm[i, k], 5.0e-3)
            ncic[i, k] = max(nc[i, k] / lcldm[i, k], 0.0)
            if settings.constant_droplet_number
                ncic[i, k] = settings.droplet_number / rho[i, k]
            end
        else
            qcic[i, k] = 0.0
            ncic[i, k] = 0.0
        end
        if qi[i, k] >= qsmall
            qiic[i, k] = min(qi[i, k] / icldm[i, k], 5.0e-3)
            niic[i, k] = max(ni[i, k] / icldm[i, k], 0.0)
            if settings.constant_ice_number
                niic[i, k] = settings.ice_number / rho[i, k]
            end
        else
            qiic[i, k] = 0.0
            niic[i, k] = 0.0
        end
        precip_frac[i, k] = cldm[i, k]
    end

    # 🌧️ precipitation fraction, top down
    if settings.precipitation_fraction_method == :in_cloud
        for i in 1:columns, k in 2:levels
            if qc[i, k] < qsmall && qi[i, k] < qsmall
                precip_frac[i, k] = precip_frac[i, k-1]
            end
        end
    else
        for i in 1:columns, k in 2:levels
            if qr[i, k-1] >= qsmall || qs[i, k-1] >= qsmall || qg[i, k-1] >= qsmall
                precip_frac[i, k] = max(precip_frac[i, k-1], precip_frac[i, k])
            end
        end
    end

    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, qcic[i, k],
            ncic[i, k], rho[i, k])
        pgam[i, k], lamc[i, k], ncic[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        qric[i, k] = qr[i, k] / precip_frac[i, k]
        nric[i, k] = nr[i, k] / precip_frac[i, k]
        qric[i, k] = min(qric[i, k], 0.01)
        if qric[i, k] < qsmall
            qric[i, k] = 0.0
            nric[i, k] = 0.0
        end
        nric[i, k] = max(nric[i, k], 0.0)
    end
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qric[i, k], nric[i, k])
        lamr[i, k], n0r[i, k], nric[i, k] = rain.slope, rain.intercept, rain.number
    end
    for k in 1:levels, i in 1:columns
        F.pgam_out[i, k] = pgam[i, k]
        F.n0r_out[i, k] = n0r[i, k]
        F.lamc_out[i, k] = lamc[i, k]
        F.lamr_out[i, k] = lamr[i, k]
    end
    return nothing
end
