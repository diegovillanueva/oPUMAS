# ⚖️🛡️ Stage 8, conservation limiters (micro_pumas_v1.F90 lines 2584–2972). Every sink of
#    a species is scaled by one ratio so the species cannot go negative in the time step;
#    then evaporation and sublimation are limited so they do not supersaturate the air.

"""
    conservation_limiters!(work, parameters)

⚖️ Stage 8. The loops and their order follow the Fortran exactly: qc (with the vapor
deposition limit), nc and rain freezing to ice, rain, ice, snow, graupel, then the
saturation limit on evaporation and sublimation.
"""
function conservation_limiters!(work, parameters)
    (; settings, saturation) = parameters
    (; cpp, rv, xxlv, xxls) = host_short_names(parameters)
    xxlv_squared, xxls_squared = xxlv^2, xxls^2
    F = work.fortran
    (; t, q, p, qvi, qc, nc, qr, nr, qi, ni, qs, ns, qg, lcldm, icldm, precip_frac) = F
    (; lamr, qcrat, dum_2D, ttmpA, qtmpAI, qvnAI, qvnA, dum1A, dum2A, dum3A) = F
    (; prc, nprc, nprc1, pra, npra, nragg, nsagg, psacws, npsacws, pracs, npracs) = F
    (; mnuccc, nnuccc, mnucct, nnucct, mnudep, nnudep, mnuccd, nnuccd, msacwi, nsacwi) = F
    (; mnuccr, nnuccr, mnuccri, nnuccri, prci, nprci, prai, nprai, pre, nsubr) = F
    (; prds, nsubs, vap_dep, ice_sublim, nsubi, vap_deps, nsubc, berg, bergs) = F
    (; pracg, npracg, psacwg, npsacwg, pgsacw, nscng, pgracs, ngracs, psacr, prdg) = F
    (; qmultg, nmultg, qmultrg, nmultrg) = F
    deltat = work.inputs.time_step
    rdeltat = 1.0 / deltat
    graupel_on = settings.do_hail || settings.do_graupel
    columns, levels = grid_size(work)

    if settings.ifs_evaporation_scaling
        for k in 1:levels, i in 1:columns
            pre[i, k] = 0.15 * pre[i, k]
        end
    end

    # 💧 conservation of qc, then the deposition limit
    for k in 1:levels, i in 1:columns
        dum = ((prc[i, k] + pra[i, k] + mnuccc[i, k] + mnucct[i, k] + msacwi[i, k] +
                psacws[i, k] + bergs[i, k] + qmultg[i, k] + psacwg[i, k] + pgsacw[i, k]) *
               lcldm[i, k] + berg[i, k]) * deltat
        if dum > qc[i, k]
            ratio = qc[i, k] * rdeltat /
                    ((prc[i, k] + pra[i, k] + mnuccc[i, k] + mnucct[i, k] + msacwi[i, k] +
                      psacws[i, k] + bergs[i, k] + qmultg[i, k] + psacwg[i, k] +
                      pgsacw[i, k]) * lcldm[i, k] + berg[i, k]) * omsm
            qmultg[i, k] = qmultg[i, k] * ratio
            psacwg[i, k] = psacwg[i, k] * ratio
            pgsacw[i, k] = pgsacw[i, k] * ratio
            prc[i, k] = prc[i, k] * ratio
            pra[i, k] = pra[i, k] * ratio
            mnuccc[i, k] = mnuccc[i, k] * ratio
            mnucct[i, k] = mnucct[i, k] * ratio
            msacwi[i, k] = msacwi[i, k] * ratio
            psacws[i, k] = psacws[i, k] * ratio
            bergs[i, k] = bergs[i, k] * ratio
            berg[i, k] = berg[i, k] * ratio
            qcrat[i, k] = ratio
        else
            qcrat[i, k] = 1.0
        end
        if qc[i, k] >= qsmall
            vap_dep[i, k] = vap_dep[i, k] * (1.0 - qcrat[i, k])
            vap_deps[i, k] = vap_deps[i, k] * (1.0 - qcrat[i, k])
        end
        mnuccd[i, k] = max(0.0, mnuccd[i, k])
        vap_dep[i, k] = max(0.0, vap_dep[i, k])
        vap_deps[i, k] = max(0.0, vap_deps[i, k])
        dum1 = vap_dep[i, k] + mnuccd[i, k] + vap_deps[i, k]
        if dum1 > 1.0e-20
            dum = (q[i, k] - qvi[i, k]) /
                  (1.0 + xxls_squared * qvi[i, k] / (cpp * rv * (t[i, k] * t[i, k]))) *
                  rdeltat
            dum = max(dum, 0.0)
            if dum1 > dum
                mnuccd[i, k] = dum * mnuccd[i, k] / dum1
                vap_dep[i, k] = dum * vap_dep[i, k] / dum1
                vap_deps[i, k] = dum * vap_deps[i, k] / dum1
            end
        end
    end

    # 💧🔢 conservation of nc; small raindrops freeze to cloud ice
    for k in 1:levels, i in 1:columns
        dum = (nprc1[i, k] + npra[i, k] + nnuccc[i, k] + nnucct[i, k] + npsacws[i, k] -
               nsubc[i, k] + npsacwg[i, k]) * lcldm[i, k] * deltat
        if dum > nc[i, k]
            ratio = nc[i, k] * rdeltat /
                    ((nprc1[i, k] + npra[i, k] + nnuccc[i, k] + nnucct[i, k] +
                      npsacws[i, k] - nsubc[i, k] + npsacwg[i, k]) * lcldm[i, k]) * omsm
            npsacwg[i, k] = npsacwg[i, k] * ratio
            nprc1[i, k] = nprc1[i, k] * ratio
            npra[i, k] = npra[i, k] * ratio
            nnuccc[i, k] = nnuccc[i, k] * ratio
            nnucct[i, k] = nnucct[i, k] * ratio
            npsacws[i, k] = npsacws[i, k] * ratio
            nsubc[i, k] = nsubc[i, k] * ratio
        end
        mnuccri[i, k] = 0.0
        nnuccri[i, k] = 0.0
        if settings.prognostic_cloud_ice && lamr[i, k] > qsmall &&
           1.0 / lamr[i, k] < settings.ice_to_snow_size
            mnuccri[i, k] = mnuccr[i, k]
            nnuccri[i, k] = nnuccr[i, k]
            mnuccr[i, k] = 0.0
            nnuccr[i, k] = 0.0
        end
    end

    # 🌧️ conservation of rain mass and number
    for k in 1:levels, i in 1:columns
        dum = ((-pre[i, k] + pracs[i, k] + mnuccr[i, k] + mnuccri[i, k] + qmultrg[i, k] +
                pracg[i, k] + pgracs[i, k]) * precip_frac[i, k] -
               (pra[i, k] + prc[i, k]) * lcldm[i, k]) * deltat
        sinks = -pre[i, k] + pracs[i, k] + mnuccr[i, k] + mnuccri[i, k] + qmultrg[i, k] +
                pracg[i, k] + pgracs[i, k]
        if dum > qr[i, k] && sinks >= qsmall
            ratio = (qr[i, k] * rdeltat + (pra[i, k] + prc[i, k]) * lcldm[i, k]) /
                    precip_frac[i, k] / sinks * omsm
            qmultrg[i, k] = qmultrg[i, k] * ratio
            pracg[i, k] = pracg[i, k] * ratio
            pgracs[i, k] = pgracs[i, k] * ratio
            pre[i, k] = pre[i, k] * ratio
            pracs[i, k] = pracs[i, k] * ratio
            mnuccr[i, k] = mnuccr[i, k] * ratio
            mnuccri[i, k] = mnuccri[i, k] * ratio
        end
        if pre[i, k] < 0.0
            nsubr[i, k] = pre[i, k] * nr[i, k] / qr[i, k]
        else
            nsubr[i, k] = 0.0
        end
        tmpnr = -nsubr[i, k] + npracs[i, k] + nnuccr[i, k] + nnuccri[i, k] - nragg[i, k] +
                npracg[i, k] + ngracs[i, k]
        dum = (tmpnr * precip_frac[i, k] - nprc[i, k] * lcldm[i, k]) * deltat
        tmpp = nr[i, k] * rdeltat + nprc[i, k] * lcldm[i, k]
        if dum > nr[i, k] && tmpnr > 0.0 && tmpp > 0.0 && precip_frac[i, k] > 0.0
            ratio = (nr[i, k] * rdeltat + nprc[i, k] * lcldm[i, k]) / precip_frac[i, k] /
                    tmpnr * omsm
            npracg[i, k] = npracg[i, k] * ratio
            ngracs[i, k] = ngracs[i, k] * ratio
            nragg[i, k] = nragg[i, k] * ratio
            npracs[i, k] = npracs[i, k] * ratio
            nnuccr[i, k] = nnuccr[i, k] * ratio
            nsubr[i, k] = nsubr[i, k] * ratio
            nnuccri[i, k] = nnuccri[i, k] * ratio
        end
    end

    # 🧊 conservation of cloud ice mass and number
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            dum = ((-mnuccc[i, k] - mnucct[i, k] - mnudep[i, k] - msacwi[i, k] -
                    qmultg[i, k]) * lcldm[i, k] + (prci[i, k] + prai[i, k]) * icldm[i, k] +
                   (-qmultrg[i, k] - mnuccri[i, k]) * precip_frac[i, k] - ice_sublim[i, k] -
                   vap_dep[i, k] - berg[i, k] - mnuccd[i, k]) * deltat
            if dum > qi[i, k]
                ratio = (qi[i, k] * rdeltat + vap_dep[i, k] + berg[i, k] + mnuccd[i, k] +
                         (mnuccc[i, k] + mnucct[i, k] + mnudep[i, k] + msacwi[i, k] +
                          qmultg[i, k]) * lcldm[i, k] +
                         (qmultrg[i, k] + mnuccri[i, k]) * precip_frac[i, k]) /
                        ((prci[i, k] + prai[i, k]) * icldm[i, k] - ice_sublim[i, k]) * omsm
                prci[i, k] = prci[i, k] * ratio
                prai[i, k] = prai[i, k] * ratio
                ice_sublim[i, k] = ice_sublim[i, k] * ratio
            end
            tmpfrz = settings.host_heterogeneous_freezing ? nnuccc[i, k] : 0.0
            dum = ((-nnucct[i, k] - tmpfrz - nnudep[i, k] - nsacwi[i, k] - nmultg[i, k]) *
                   lcldm[i, k] + (nprci[i, k] + nprai[i, k] - nsubi[i, k]) * icldm[i, k] +
                   (-nmultrg[i, k] - nnuccri[i, k]) * precip_frac[i, k] - nnuccd[i, k]) *
                  deltat
            if dum > ni[i, k]
                ratio = (ni[i, k] * rdeltat + nnuccd[i, k] +
                         (nnucct[i, k] + tmpfrz + nnudep[i, k] + nsacwi[i, k] +
                          nmultg[i, k]) * lcldm[i, k] +
                         (nnuccri[i, k] + nmultrg[i, k]) * precip_frac[i, k]) /
                        ((nprci[i, k] + nprai[i, k] - nsubi[i, k]) * icldm[i, k]) * omsm
                nprci[i, k] = nprci[i, k] * ratio
                nprai[i, k] = nprai[i, k] * ratio
                nsubi[i, k] = nsubi[i, k] * ratio
            end
        end
    end

    # ❄️ conservation of snow mass and number
    for k in 1:levels, i in 1:columns
        if graupel_on
            dum = (-(prds[i, k] + pracs[i, k] - psacr[i, k]) * precip_frac[i, k] -
                   (prai[i, k] + prci[i, k]) * icldm[i, k] -
                   (bergs[i, k] + psacws[i, k]) * lcldm[i, k] - vap_deps[i, k]) * deltat
        else
            dum = (-(prds[i, k] + pracs[i, k] + mnuccr[i, k]) * precip_frac[i, k] -
                   (prai[i, k] + prci[i, k]) * icldm[i, k] -
                   (bergs[i, k] + psacws[i, k]) * lcldm[i, k] - vap_deps[i, k]) * deltat
        end
        if dum > qs[i, k] && (psacr[i, k] - prds[i, k]) >= qsmall
            if graupel_on
                ratio = (qs[i, k] * rdeltat + (prai[i, k] + prci[i, k]) * icldm[i, k] +
                         (bergs[i, k] + psacws[i, k]) * lcldm[i, k] + vap_deps[i, k] +
                         pracs[i, k] * precip_frac[i, k]) /
                        precip_frac[i, k] / (psacr[i, k] - prds[i, k]) * omsm
                psacr[i, k] = psacr[i, k] * ratio
            else
                ratio = (qs[i, k] * rdeltat + (prai[i, k] + prci[i, k]) * icldm[i, k] +
                         (bergs[i, k] + psacws[i, k]) * lcldm[i, k] + vap_deps[i, k] +
                         (pracs[i, k] + mnuccr[i, k]) * precip_frac[i, k]) /
                        precip_frac[i, k] / (-prds[i, k]) * omsm
            end
            prds[i, k] = prds[i, k] * ratio
        end
        nsubs[i, k] = 0.0
        if graupel_on
            dum = ((-nsagg[i, k] - nsubs[i, k] + ngracs[i, k]) * precip_frac[i, k] -
                   nprci[i, k] * icldm[i, k] + nscng[i, k] * lcldm[i, k]) * deltat
        else
            dum = ((-nsagg[i, k] - nsubs[i, k] - nnuccr[i, k]) * precip_frac[i, k] -
                   nprci[i, k] * icldm[i, k]) * deltat
        end
        if dum > ns[i, k]
            if graupel_on
                ratio = (ns[i, k] * rdeltat + nprci[i, k] * icldm[i, k]) /
                        precip_frac[i, k] /
                        (-nsubs[i, k] - nsagg[i, k] + ngracs[i, k] +
                         lcldm[i, k] / precip_frac[i, k] * nscng[i, k]) * omsm
                nscng[i, k] = nscng[i, k] * ratio
                ngracs[i, k] = ngracs[i, k] * ratio
            else
                ratio = (ns[i, k] * rdeltat + nnuccr[i, k] * precip_frac[i, k] +
                         nprci[i, k] * icldm[i, k]) / precip_frac[i, k] /
                        (-nsubs[i, k] - nsagg[i, k]) * omsm
            end
            nsubs[i, k] = nsubs[i, k] * ratio
            nsagg[i, k] = nsagg[i, k] * ratio
        end
    end

    # 🟤 conservation of graupel mass (graupel number has no sinks)
    if graupel_on
        for k in 1:levels, i in 1:columns
            dum = ((-pracg[i, k] - pgracs[i, k] - prdg[i, k] - psacr[i, k] - mnuccr[i, k]) *
                   precip_frac[i, k] + (-psacwg[i, k] - pgsacw[i, k]) * lcldm[i, k]) * deltat
            if dum > qg[i, k]
                ratio = (qg[i, k] * rdeltat +
                         (pracg[i, k] + pgracs[i, k] + psacr[i, k] + mnuccr[i, k]) *
                         precip_frac[i, k] + (psacwg[i, k] + pgsacw[i, k]) * lcldm[i, k]) /
                        ((-prdg[i, k]) * precip_frac[i, k]) * omsm
                prdg[i, k] = prdg[i, k] * ratio
            end
        end
    end

    # 🌫️ limit evaporation and sublimation so that they do not supersaturate the air
    for k in 1:levels, i in 1:columns
        qtmpAI[i, k] = q[i, k] - (ice_sublim[i, k] + vap_dep[i, k] + mnuccd[i, k] +
                                  vap_deps[i, k] +
                                  (pre[i, k] + prds[i, k] + prdg[i, k]) * precip_frac[i, k]) *
                                 deltat
        ttmpA[i, k] = t[i, k] + ((pre[i, k] * precip_frac[i, k]) * xxlv +
                                 ((prds[i, k] + prdg[i, k]) * precip_frac[i, k] +
                                  vap_dep[i, k] + vap_deps[i, k] + ice_sublim[i, k] +
                                  mnuccd[i, k]) * xxls) * deltat / cpp
    end
    for k in 1:levels, i in 1:columns
        qvnAI[i, k] = saturation_over_water(saturation, ttmpA[i, k], p[i, k]).specific_humidity
    end
    evaporating(i, k) =
        (pre[i, k] + prds[i, k] + prdg[i, k]) * precip_frac[i, k] + ice_sublim[i, k] <
        -1.0e-20 && qtmpAI[i, k] > qvnAI[i, k]
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            total = (pre[i, k] + prds[i, k] + prdg[i, k]) * precip_frac[i, k] +
                    ice_sublim[i, k]
            dum1A[i, k] = pre[i, k] * precip_frac[i, k] / total
            dum2A[i, k] = prds[i, k] * precip_frac[i, k] / total
            dum3A[i, k] = prdg[i, k] * precip_frac[i, k] / total
            ttmpA[i, k] = t[i, k] +
                          ((vap_dep[i, k] + vap_deps[i, k] + mnuccd[i, k]) * xxls) *
                          deltat / cpp
            dum_2D[i, k] = q[i, k] - (vap_dep[i, k] + vap_deps[i, k] + mnuccd[i, k]) * deltat
        end
    end
    for k in 1:levels, i in 1:columns
        qvnA[i, k] = saturation_over_water(saturation, ttmpA[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            dum = (dum_2D[i, k] - qvnA[i, k]) /
                  (1.0 + xxlv_squared * qvnA[i, k] / (cpp * rv * (ttmpA[i, k] * ttmpA[i, k])))
            dum = min(dum, 0.0)
            pre[i, k] = dum * dum1A[i, k] * rdeltat / precip_frac[i, k]
        end
    end
    for k in 1:levels, i in 1:columns
        qvnA[i, k] = saturation_over_ice(saturation, ttmpA[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            dum = (dum_2D[i, k] - qvnA[i, k]) /
                  (1.0 + xxls_squared * qvnA[i, k] / (cpp * rv * (ttmpA[i, k] * ttmpA[i, k])))
            dum = min(dum, 0.0)
            prds[i, k] = dum * dum2A[i, k] * rdeltat / precip_frac[i, k]
            prdg[i, k] = dum * dum3A[i, k] * rdeltat / precip_frac[i, k]
            dum1A[i, k] = (1.0 - dum1A[i, k] - dum2A[i, k] - dum3A[i, k])
            ice_sublim[i, k] = dum * dum1A[i, k] * rdeltat
        end
    end
    return nothing
end
