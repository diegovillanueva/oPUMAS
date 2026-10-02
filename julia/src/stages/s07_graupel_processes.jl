# 🟤🌧️ Stage 7, process rates of graupel or hail and evaporation of precipitation
#    (micro_pumas_v1.F90 lines 2531–2583).

"""
    graupel_processes!(work, parameters)

🟤 Stage 7. Collection of snow, cloud liquid and rain by graupel, rime splintering, and
evaporation and sublimation of rain, snow and graupel.
"""
function graupel_processes!(work, parameters)
    (; settings) = parameters
    if !(settings.do_hail || settings.do_graupel)
        error("precipitation evaporation without graupel or hail is not translated yet")
    end
    dense = dense_precipitation(parameters)
    F = work.fortran
    (; t, rho, arn, asn, agn, lcldm, precip_frac) = F
    (; qcic, ncic, qiic, qric, nric, lamr, n0r, unr, umr) = F
    (; qsic, nsic, lams, n0s, ums, qgic, lamg, n0g, ung, umg) = F
    (; psacr, psacwg, npsacwg, psacws, pgsacw, nscng, pracg, npracg, pracs, npracs) = F
    (; pgracs, ngracs, qmultg, nmultg, qmultrg, nmultrg, pre, prds, prdg, am_evp_st) = F
    deltat = work.inputs.time_step
    iaccr = settings.ice_accretion_factor
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        psacr[i, k] = snow_collection_by_rain(umr[i, k], ums[i, k], rho[i, k], qsic[i, k],
            qric[i, k], lamr[i, k], n0r[i, k], lams[i, k], n0s[i, k])
    end
    for k in 1:levels, i in 1:columns
        rate = droplet_collection_by_graupel(qgic[i, k], qcic[i, k], ncic[i, k], rho[i, k],
            n0g[i, k], lamg[i, k], dense.b, agn[i, k])
        psacwg[i, k], npsacwg[i, k] = rate.mass, rate.number
    end
    for k in 1:levels, i in 1:columns
        psacwg[i, k] = psacwg[i, k] * iaccr
        npsacwg[i, k] = npsacwg[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        rate = rimed_snow_to_graupel(parameters, psacws[i, k], qsic[i, k], qcic[i, k],
            nsic[i, k], rho[i, k], rhosn, dense.density, asn[i, k], lams[i, k], n0s[i, k],
            deltat)
        pgsacw[i, k], nscng[i, k], psacws[i, k] = rate.mass, rate.number, rate.riming
    end
    for k in 1:levels, i in 1:columns
        rate = rain_collection_by_graupel(qric[i, k], qgic[i, k], umg[i, k], umr[i, k],
            ung[i, k], unr[i, k], rho[i, k], n0r[i, k], lamr[i, k], n0g[i, k], lamg[i, k])
        pracg[i, k], npracg[i, k] = rate.mass, rate.number
    end
    for k in 1:levels, i in 1:columns
        pracg[i, k] = pracg[i, k] * iaccr
        npracg[i, k] = npracg[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        rate = rain_snow_collisions_to_graupel(pracs[i, k], npracs[i, k], psacr[i, k],
            qsic[i, k], qric[i, k], nric[i, k], nsic[i, k], n0s[i, k], lams[i, k],
            n0r[i, k], lamr[i, k], deltat)
        pgracs[i, k], ngracs[i, k] = rate.mass, rate.number
        pracs[i, k], npracs[i, k] = rate.rain_by_snow, rate.rain_number_by_snow
        psacr[i, k] = rate.snow_by_rain
    end
    for k in 1:levels, i in 1:columns
        rate = rime_splintering_by_graupel(parameters, t[i, k], qcic[i, k], qric[i, k],
            qgic[i, k], psacwg[i, k], pracg[i, k])
        qmultg[i, k], nmultg[i, k] = rate.droplet_mass, rate.droplet_number
        qmultrg[i, k], nmultrg[i, k] = rate.rain_mass, rate.rain_number
        psacwg[i, k], pracg[i, k] = rate.droplet_riming, rate.rain_riming
    end
    for k in 1:levels, i in 1:columns
        rain = (q = qric[i, k], λ = lamr[i, k], N₀ = n0r[i, k], a = arn[i, k])
        snow = (q = qsic[i, k], λ = lams[i, k], N₀ = n0s[i, k], a = asn[i, k])
        graupel = (q = qgic[i, k], λ = lamg[i, k], N₀ = n0g[i, k], a = agn[i, k],
            b = dense.b)
        rate = precipitation_evaporation(parameters, air_at(F, i, k), lcldm[i, k],
            precip_frac[i, k], qcic[i, k], qiic[i, k], rain, snow, graupel)
        pre[i, k], prds[i, k], prdg[i, k] = rate.rain, rate.snow, rate.graupel
        am_evp_st[i, k] = rate.area
    end
    return nothing
end
