# ❄️🧊 Stage 6, process rates of ice, snow and mixed phase (micro_pumas_v1.F90 lines
#    2166–2530).

"Scalar air NamedTuple of the vapor diffusion twins at point (i, k)."
function air_at(F, i, k)
    return (T = F.t[i, k], ρ = F.rho[i, k], diffusivity = F.dv[i, k],
        viscosity = F.mu[i, k], schmidt = F.sc[i, k], q_vapor = F.q[i, k],
        q_sat_liquid = F.qvl[i, k], q_sat_ice = F.qvi[i, k])
end

"""
    ice_snow_processes!(work, parameters)

❄️ Stage 6. Size distributions of ice, snow and graupel with provisional fall speeds,
then ice autoconversion, droplet freezing, aggregation, riming, rain freezing, accretion,
the Bergeron process and vapor deposition.
"""
function ice_snow_processes!(work, parameters)
    (; settings, hydrometeors, gamma) = parameters
    dense = dense_precipitation(parameters)
    F = work.fortran
    (; t, q, rho, mu, rhof, arn, asn, agn, relvar, accre_enhan, cldm, lcldm, icldm) = F
    (; qi, ni, qs, ns, qg, ng, precip_frac, dum_2D) = F
    (; qcic, ncic, pgam, lamc, qiic, niic, lami, n0i, qric, nric, lamr, n0r, unr) = F
    (; qsic, nsic, lams, n0s, uns, qgic, ngic, lamg, n0g, ung, umr, ums, umg) = F
    (; prci, nprci, prc, nprc, pra, npra, nragg, nsagg, psacws, npsacws, msacwi, nsacwi) = F
    (; pracs, npracs, mnuccr, nnuccr, mnuccc, nnuccc, mnucct, nnucct, mnudep, nnudep) = F
    (; prai, nprai, bergs, vap_deps, vap_dep, ice_sublim, berg, nsubi, nsubc) = F
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        ice = size_distribution(hydrometeors.cloud_ice, qiic[i, k], niic[i, k])
        lami[i, k], n0i[i, k], niic[i, k] = ice.slope, ice.intercept, ice.number
    end

    # 🧊➡️❄️ autoconversion of cloud ice to snow
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            rate = ice_to_snow_autoconversion(parameters, t[i, k], qiic[i, k], lami[i, k],
                n0i[i, k])
            prci[i, k], nprci[i, k] = rate.mass, rate.number
        else
            prci[i, k] = F.tnd_qsnow[i, k] / cldm[i, k]
            nprci[i, k] = F.tnd_nsnow[i, k] / cldm[i, k]
        end
    end

    # ❄️🟤 in-precipitation snow and graupel
    for k in 1:levels, i in 1:columns
        qsic[i, k] = qs[i, k] / precip_frac[i, k]
        nsic[i, k] = ns[i, k] / precip_frac[i, k]
        qsic[i, k] = min(qsic[i, k], 0.01)
        if qsic[i, k] < qsmall
            qsic[i, k] = 0.0
            nsic[i, k] = 0.0
        end
        nsic[i, k] = max(nsic[i, k], 0.0)
        qgic[i, k] = qg[i, k] / precip_frac[i, k]
        ngic[i, k] = ng[i, k] / precip_frac[i, k]
        qgic[i, k] = min(qgic[i, k], 0.01)
        if qgic[i, k] < qsmall
            qgic[i, k] = 0.0
            ngic[i, k] = 0.0
        end
        ngic[i, k] = max(ngic[i, k], 0.0)
    end

    # 🌧️ provisional rain fall speeds (m/s)
    br = fall_speed_law.rain.b
    for k in 1:levels, i in 1:columns
        if lamr[i, k] >= qsmall
            dum_2D[i, k] = lamr[i, k]^br
            unr[i, k] = min(arn[i, k] * gamma.rain_plus1 / dum_2D[i, k], 9.1 * rhof[i, k])
            umr[i, k] = min(arn[i, k] * gamma.rain_plus4 / (6.0 * dum_2D[i, k]),
                9.1 * rhof[i, k])
        else
            umr[i, k] = 0.0
            unr[i, k] = 0.0
        end
    end

    # ❄️ snow size distribution and provisional fall speeds
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, qsic[i, k], nsic[i, k])
        lams[i, k], n0s[i, k], nsic[i, k] = snow.slope, snow.intercept, snow.number
    end
    bs = fall_speed_law.snow.b
    for k in 1:levels, i in 1:columns
        if settings.ifs_snow_fall_speed
            ums[i, k] = lams[i, k] > 0.0 ? 1.0 : 0.0
            uns[i, k] = lams[i, k] > 0.0 ? 1.0 : 0.0
        elseif lams[i, k] > 0.0
            dum_2D[i, k] = lams[i, k]^bs
            ums[i, k] = min(asn[i, k] * gamma.snow_plus4 / (6.0 * dum_2D[i, k]),
                1.2 * rhof[i, k])
            ums[i, k] = ums[i, k] * settings.snow_fall_speed_factor
            uns[i, k] = min(asn[i, k] * gamma.snow_plus1 / dum_2D[i, k], 1.2 * rhof[i, k])
        else
            ums[i, k] = 0.0
            uns[i, k] = 0.0
        end
    end

    # 🟤 graupel (or hail) size distribution and provisional fall speeds
    if settings.do_hail || settings.do_graupel
        dense_props = settings.do_hail ? hydrometeors.hail : hydrometeors.graupel
        for k in 1:levels, i in 1:columns
            graupel = size_distribution(dense_props, qgic[i, k], ngic[i, k])
            lamg[i, k], n0g[i, k], ngic[i, k] =
                graupel.slope, graupel.intercept, graupel.number
        end
    end
    for k in 1:levels, i in 1:columns
        if lamg[i, k] > 0.0
            dum_2D[i, k] = lamg[i, k]^dense.b
            umg[i, k] = min(agn[i, k] * gamma.graupel_plus4 / (6.0 * dum_2D[i, k]),
                20.0 * rhof[i, k])
            ung[i, k] = min(agn[i, k] * gamma.graupel_plus1 / dum_2D[i, k],
                20.0 * rhof[i, k])
        else
            umg[i, k] = 0.0
            ung[i, k] = 0.0
        end
    end

    # 💧➡️🧊 droplet freezing from the host's heterogeneous freezing rates
    if settings.prognostic_cloud_ice
        settings.host_heterogeneous_freezing ||
            error("Bigg droplet freezing (hetfrz off) is not translated yet")
        for k in 1:levels, i in 1:columns
            mi0l = qcic[i, k] / max(ncic[i, k], 1.0e6 / rho[i, k])
            mi0l = max(mi0l_min, mi0l)
            if qcic[i, k] >= qsmall
                nnuccc[i, k] = F.frzimm[i, k] * 1.0e6 / rho[i, k]
                mnuccc[i, k] = nnuccc[i, k] * mi0l
                nnucct[i, k] = F.frzcnt[i, k] * 1.0e6 / rho[i, k]
                mnucct[i, k] = nnucct[i, k] * mi0l
                nnudep[i, k] = F.frzdep[i, k] * 1.0e6 / rho[i, k]
                mnudep[i, k] = nnudep[i, k] * mi0
            else
                nnuccc[i, k] = 0.0
                mnuccc[i, k] = 0.0
                nnucct[i, k] = 0.0
                mnucct[i, k] = 0.0
                nnudep[i, k] = 0.0
                mnudep[i, k] = 0.0
            end
        end
    else
        for array in (mnuccc, nnuccc, mnucct, nnucct, mnudep, nnudep)
            fill!(array, 0.0)
        end
    end

    iaccr = settings.ice_accretion_factor
    for k in 1:levels, i in 1:columns
        nsagg[i, k] = snow_self_aggregation(parameters, t[i, k], rho[i, k], asn[i, k],
            rhosn, qsic[i, k], nsic[i, k])
    end
    for k in 1:levels, i in 1:columns
        riming = droplet_collection_by_snow(parameters, t[i, k], rho[i, k], asn[i, k],
            uns[i, k], mu[i, k], qcic[i, k], ncic[i, k], qsic[i, k], pgam[i, k],
            lamc[i, k], lams[i, k], n0s[i, k])
        psacws[i, k], npsacws[i, k] = riming.mass, riming.number
    end
    for k in 1:levels, i in 1:columns
        psacws[i, k] = psacws[i, k] * iaccr
        npsacws[i, k] = npsacws[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            splinters = rime_splintering_by_snow(t[i, k], psacws[i, k])
            msacwi[i, k], nsacwi[i, k] = splinters.mass, splinters.number
            psacws[i, k] = splinters.riming
        else
            nsacwi[i, k] = 0.0
            msacwi[i, k] = 0.0
        end
    end
    for k in 1:levels, i in 1:columns
        collection = rain_collection_by_snow(parameters, t[i, k], rho[i, k], umr[i, k],
            ums[i, k], unr[i, k], uns[i, k], qric[i, k], qsic[i, k], lamr[i, k], n0r[i, k],
            lams[i, k], n0s[i, k])
        pracs[i, k], npracs[i, k] = collection.mass, collection.number
    end
    for k in 1:levels, i in 1:columns
        freezing = rain_immersion_freezing(parameters, t[i, k], qric[i, k], nric[i, k],
            lamr[i, k])
        mnuccr[i, k], nnuccr[i, k] = freezing.mass, freezing.number
    end

    # 💧🌧️ accretion of cloud liquid by rain
    if settings.warm_rain == :sb2001
        (; qctend_SB2001, nctend_SB2001, qrtend_SB2001, nrtend_SB2001) = F
        for k in 1:levels, i in 1:columns
            rate = sb2001_accretion(qcic[i, k], ncic[i, k], qric[i, k], rho[i, k])
            pra[i, k], npra[i, k] = rate.mass, rate.droplet_number
        end
        for k in 1:levels, i in 1:columns
            nctend_SB2001[i, k] = nctend_SB2001[i, k] + npra[i, k]
            qctend_SB2001[i, k] = qctend_SB2001[i, k] + pra[i, k]
            nrtend_SB2001[i, k] = nrtend_SB2001[i, k] + npra[i, k]
            qrtend_SB2001[i, k] = qrtend_SB2001[i, k] - pra[i, k]
        end
    end
    if settings.warm_rain == :kk2000
        (; qctend_KK2000, nctend_KK2000, qrtend_KK2000, nrtend_KK2000) = F
        for k in 1:levels, i in 1:columns
            rtmp = qric[i, k]
            ctmp = qcic[i, k]
            ntmp = ncic[i, k]
            if settings.accretion_sees_autoconversion
                rtmp = rtmp + prc[i, k] * deltat
                ctmp = ctmp - prc[i, k] * deltat
                ntmp = ntmp - nprc[i, k] * deltat
            end
            rate = accretion(parameters, rtmp, ctmp, ntmp, relvar[i, k], accre_enhan[i, k])
            pra[i, k], npra[i, k] = rate.mass, rate.droplet_number
        end
        for k in 1:levels, i in 1:columns
            nctend_KK2000[i, k] = nctend_KK2000[i, k] + npra[i, k]
            qctend_KK2000[i, k] = qctend_KK2000[i, k] + pra[i, k]
            nrtend_KK2000[i, k] = nrtend_KK2000[i, k] + npra[i, k]
            qrtend_KK2000[i, k] = qrtend_KK2000[i, k] - pra[i, k]
        end
    end
    for k in 1:levels, i in 1:columns
        pra[i, k] = pra[i, k] * settings.accretion_enhancement
        npra[i, k] = npra[i, k] * settings.accretion_enhancement
    end
    for k in 1:levels, i in 1:columns
        nragg[i, k] = rain_self_collection(rho[i, k], qric[i, k], nric[i, k])
    end

    # 🧊❄️ accretion of cloud ice by snow
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            rate = ice_collection_by_snow(parameters, t[i, k], rho[i, k], asn[i, k],
                qiic[i, k], niic[i, k], qsic[i, k], lams[i, k], n0s[i, k])
            prai[i, k], nprai[i, k] = rate.mass, rate.number
        else
            prai[i, k] = 0.0
            nprai[i, k] = 0.0
        end
    end

    # 🌫️ Bergeron process and vapor deposition onto snow and ice
    berg_eff = settings.bergeron_efficiency
    for k in 1:levels, i in 1:columns
        bergs[i, k] = bergeron_to_snow(parameters, air_at(F, i, k), asn[i, k], qcic[i, k],
            qsic[i, k], lams[i, k], n0s[i, k])
    end
    for k in 1:levels, i in 1:columns
        bergs[i, k] = bergs[i, k] * berg_eff
    end
    for k in 1:levels, i in 1:columns
        vap_deps[i, k] = snow_deposition(parameters, air_at(F, i, k), qs[i, k], ns[i, k],
            precip_frac[i, k], asn[i, k])
    end
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            rate = ice_deposition(parameters, air_at(F, i, k), qi[i, k], ni[i, k],
                icldm[i, k])
            berg[i, k], vap_dep[i, k] = rate.bergeron, rate.deposition
            ice_sublim[i, k] = rate.sublimation
        end
        for k in 1:levels, i in 1:columns
            berg[i, k] = berg[i, k] * berg_eff
            if ice_sublim[i, k] < 0.0 && qi[i, k] > qsmall && icldm[i, k] > mincld
                nsubi[i, k] = sublim_factor * ice_sublim[i, k] / qi[i, k] * ni[i, k] /
                              icldm[i, k]
            else
                nsubi[i, k] = 0.0
            end
            nsubc[i, k] = 0.0
        end
    end
    return nothing
end
