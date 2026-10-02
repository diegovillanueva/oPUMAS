# 📈🧮 Stage 9, tendencies of the conversion processes (micro_pumas_v1.F90 lines
#    2972–3172): local rates times their cloud or precipitation fraction become grid-mean
#    tendencies, history process rates and host diagnostics.

"""
    process_tendencies!(work, parameters)

📈 Stage 9. Adds the limited process rates to the tendencies (which already hold the
instantaneous melting and freezing of stage 3), fills the grid-mean process rates and
diagnostics, and recomputes the rain distribution and `rercld`.
"""
function process_tendencies!(work, parameters)
    (; settings, hydrometeors) = parameters
    (; xxlv, xxls, xlf) = host_short_names(parameters)
    F = work.fortran
    (; qc, nc, qr, nr, qs, ns, qg, ng, rho, lcldm, icldm, precip_frac) = F
    (; qvlat, tlat, qctend, nctend, qitend, nitend, qrtend, nrtend, qstend, nstend) = F
    (; qgtend, ngtend) = F
    (; prc, nprc, nprc1, pra, npra, nragg, nsagg, psacws, npsacws, pracs, npracs) = F
    (; mnuccc, nnuccc, mnucct, nnucct, mnudep, nnudep, mnuccd, nnuccd, msacwi, nsacwi) = F
    (; mnuccr, nnuccr, mnuccri, nnuccri, prci, nprci, prai, nprai, pre, nsubr) = F
    (; prds, nsubs, vap_dep, ice_sublim, nsubi, vap_deps, nsubc, berg, bergs) = F
    (; pracg, npracg, psacwg, npsacwg, pgsacw, nscng, pgracs, ngracs, psacr, prdg) = F
    (; qmultg, nmultg, qmultrg, nmultrg, ninstsm, ninstgm) = F
    (; cmeout, nevapr, prer_evap, prain, qcsinksum_rate1ord) = F
    (; qrout, nrout, qsout, nsout, qgout, ngout) = F
    (; qric, nric, lamr, n0r, lamc, pgam, qcic, ncic, rercld) = F
    deltat = work.inputs.time_step
    graupel_on = settings.do_hail || settings.do_graupel
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        qvlat[i, k] = qvlat[i, k] - (pre[i, k] + prds[i, k]) * precip_frac[i, k] -
                      vap_dep[i, k] - vap_deps[i, k] - ice_sublim[i, k] - mnuccd[i, k] -
                      mnudep[i, k] * lcldm[i, k] - prdg[i, k] * precip_frac[i, k]
        tlat[i, k] = tlat[i, k] +
                     ((pre[i, k] * precip_frac[i, k]) * xxlv +
                      ((prds[i, k] + prdg[i, k]) * precip_frac[i, k] + vap_dep[i, k] +
                       vap_deps[i, k] + ice_sublim[i, k] + mnuccd[i, k] +
                       mnudep[i, k] * lcldm[i, k]) * xxls +
                      ((bergs[i, k] + psacws[i, k] + mnuccc[i, k] + mnucct[i, k] +
                        msacwi[i, k] + psacwg[i, k] + qmultg[i, k] + pgsacw[i, k]) *
                       lcldm[i, k] +
                       (mnuccr[i, k] + pracs[i, k] + mnuccri[i, k] + pracg[i, k] +
                        pgracs[i, k] + qmultrg[i, k]) * precip_frac[i, k] +
                       berg[i, k]) * xlf)
        qctend[i, k] = qctend[i, k] +
                       (-pra[i, k] - prc[i, k] - mnuccc[i, k] - mnucct[i, k] -
                        msacwi[i, k] - psacws[i, k] - bergs[i, k] - qmultg[i, k] -
                        psacwg[i, k] - pgsacw[i, k]) * lcldm[i, k] - berg[i, k]
        if settings.prognostic_cloud_ice
            qitend[i, k] = qitend[i, k] +
                           (mnuccc[i, k] + mnucct[i, k] + mnudep[i, k] + msacwi[i, k] +
                            qmultg[i, k]) * lcldm[i, k] +
                           (-prci[i, k] - prai[i, k]) * icldm[i, k] + vap_dep[i, k] +
                           berg[i, k] + ice_sublim[i, k] + mnuccd[i, k] +
                           (mnuccri[i, k] + qmultrg[i, k]) * precip_frac[i, k]
        end
        qrtend[i, k] = qrtend[i, k] + (pra[i, k] + prc[i, k]) * lcldm[i, k] +
                       (pre[i, k] - pracs[i, k] - mnuccr[i, k] - mnuccri[i, k] -
                        qmultrg[i, k] - pracg[i, k] - pgracs[i, k]) * precip_frac[i, k]
        if graupel_on
            qgtend[i, k] = qgtend[i, k] +
                           (pracg[i, k] + pgracs[i, k] + prdg[i, k] + psacr[i, k] +
                            mnuccr[i, k]) * precip_frac[i, k] +
                           (psacwg[i, k] + pgsacw[i, k]) * lcldm[i, k]
            qstend[i, k] = qstend[i, k] + (prai[i, k] + prci[i, k]) * icldm[i, k] +
                           (psacws[i, k] + bergs[i, k]) * lcldm[i, k] +
                           (prds[i, k] + pracs[i, k] - psacr[i, k]) * precip_frac[i, k] +
                           vap_deps[i, k]
        else
            qstend[i, k] = qstend[i, k] + (prai[i, k] + prci[i, k]) * icldm[i, k] +
                           (psacws[i, k] + bergs[i, k]) * lcldm[i, k] +
                           (prds[i, k] + pracs[i, k] + mnuccr[i, k]) * precip_frac[i, k] +
                           vap_deps[i, k]
        end
    end

    for k in 1:levels, i in 1:columns
        cmeout[i, k] = vap_dep[i, k] + ice_sublim[i, k] + mnuccd[i, k] + vap_deps[i, k]
        F.cmeitot[i, k] = vap_dep[i, k] + ice_sublim[i, k] + mnuccd[i, k] + vap_deps[i, k]
        F.evapsnow[i, k] = (-prds[i, k] - prdg[i, k]) * precip_frac[i, k]
        nevapr[i, k] = -pre[i, k] * precip_frac[i, k]
        prer_evap[i, k] = -pre[i, k] * precip_frac[i, k]
        prain[i, k] = (pra[i, k] + prc[i, k]) * lcldm[i, k] +
                      (-pracs[i, k] - mnuccr[i, k] - mnuccri[i, k]) * precip_frac[i, k]
        if graupel_on
            F.prodsnow[i, k] = (prai[i, k] + prci[i, k]) * icldm[i, k] +
                               (psacws[i, k] + bergs[i, k]) * lcldm[i, k] +
                               (pracs[i, k]) * precip_frac[i, k] + vap_deps[i, k]
        else
            F.prodsnow[i, k] = (prai[i, k] + prci[i, k]) * icldm[i, k] +
                               (psacws[i, k] + bergs[i, k]) * lcldm[i, k] +
                               (pracs[i, k] + mnuccr[i, k]) * precip_frac[i, k] +
                               vap_deps[i, k]
        end
        qcsinksum_rate1ord[i, k] = (pra[i, k] + prc[i, k] + psacws[i, k] + psacwg[i, k] +
                                    pgsacw[i, k]) * lcldm[i, k]
        qcsinksum_rate1ord[i, k] = qcsinksum_rate1ord[i, k] / max(qc[i, k], 1.0e-30)
    end

    for k in 1:levels, i in 1:columns
        F.pratot[i, k] = pra[i, k] * lcldm[i, k]
        F.prctot[i, k] = prc[i, k] * lcldm[i, k]
        F.mnuccctot[i, k] = mnuccc[i, k] * lcldm[i, k]
        F.mnudeptot[i, k] = mnudep[i, k] * lcldm[i, k]
        F.mnuccttot[i, k] = mnucct[i, k] * lcldm[i, k]
        F.msacwitot[i, k] = msacwi[i, k] * lcldm[i, k]
        F.psacwstot[i, k] = psacws[i, k] * lcldm[i, k]
        F.bergstot[i, k] = bergs[i, k] * lcldm[i, k]
        F.vapdepstot[i, k] = vap_deps[i, k]
        F.bergtot[i, k] = berg[i, k]
        F.prcitot[i, k] = prci[i, k] * icldm[i, k]
        F.praitot[i, k] = prai[i, k] * icldm[i, k]
        F.mnuccdtot[i, k] = mnuccd[i, k] * icldm[i, k]
        F.pracstot[i, k] = pracs[i, k] * precip_frac[i, k]
        F.mnuccrtot[i, k] = mnuccr[i, k] * precip_frac[i, k]
        F.mnuccritot[i, k] = mnuccri[i, k] * precip_frac[i, k]
        F.psacrtot[i, k] = psacr[i, k] * precip_frac[i, k]
        F.pracgtot[i, k] = pracg[i, k] * precip_frac[i, k]
        F.psacwgtot[i, k] = psacwg[i, k] * lcldm[i, k]
        F.pgsacwtot[i, k] = pgsacw[i, k] * lcldm[i, k]
        F.pgracstot[i, k] = pgracs[i, k] * precip_frac[i, k]
        F.prdgtot[i, k] = prdg[i, k] * precip_frac[i, k]
        F.qmultgtot[i, k] = qmultg[i, k] * lcldm[i, k]
        F.qmultrgtot[i, k] = qmultrg[i, k] * precip_frac[i, k]
        F.npracgtot[i, k] = npracg[i, k] * precip_frac[i, k]
        F.nscngtot[i, k] = nscng[i, k] * lcldm[i, k]
        F.ngracstot[i, k] = ngracs[i, k] * precip_frac[i, k]
        F.nmultgtot[i, k] = nmultg[i, k] * lcldm[i, k]
        F.nmultrgtot[i, k] = nmultrg[i, k] * precip_frac[i, k]
        F.npsacwgtot[i, k] = npsacwg[i, k] * lcldm[i, k]
        F.nnuccctot[i, k] = nnuccc[i, k] * lcldm[i, k]
        F.nnuccttot[i, k] = nnucct[i, k] * lcldm[i, k]
        F.nnuccdtot[i, k] = nnuccd[i, k] * icldm[i, k]
        F.nnudeptot[i, k] = nnudep[i, k] * lcldm[i, k]
        F.nnuccrtot[i, k] = nnuccr[i, k] * precip_frac[i, k]
        F.nnuccritot[i, k] = nnuccri[i, k] * precip_frac[i, k]
        F.nsacwitot[i, k] = nsacwi[i, k] * lcldm[i, k]
        F.npratot[i, k] = npra[i, k] * lcldm[i, k]
        F.npsacwstot[i, k] = npsacws[i, k] * lcldm[i, k]
        F.npraitot[i, k] = nprai[i, k] * icldm[i, k]
        F.npracstot[i, k] = npracs[i, k] * precip_frac[i, k]
        F.nprctot[i, k] = nprc[i, k] * lcldm[i, k]
        F.nraggtot[i, k] = nragg[i, k] * precip_frac[i, k]
        F.nprcitot[i, k] = nprci[i, k] * icldm[i, k]
        F.nmeltstot[i, k] = ninstsm[i, k] / deltat
        F.nmeltgtot[i, k] = ninstgm[i, k] / deltat
    end

    for k in 1:levels, i in 1:columns
        nctend[i, k] = nctend[i, k] +
                       (-nnuccc[i, k] - nnucct[i, k] - npsacws[i, k] + nsubc[i, k] -
                        npra[i, k] - nprc1[i, k] - npsacwg[i, k]) * lcldm[i, k]
        if settings.prognostic_cloud_ice
            tmpfrz = settings.host_heterogeneous_freezing ? nnuccc[i, k] : 0.0
            nitend[i, k] = nitend[i, k] + nnuccd[i, k] +
                           (nnucct[i, k] + tmpfrz + nnudep[i, k] + nsacwi[i, k] +
                            nmultg[i, k]) * lcldm[i, k] +
                           (nsubi[i, k] - nprci[i, k] - nprai[i, k]) * icldm[i, k] +
                           (nnuccri[i, k] + nmultrg[i, k]) * precip_frac[i, k]
        end
        if graupel_on
            nstend[i, k] = nstend[i, k] +
                           (nsubs[i, k] + nsagg[i, k] - ngracs[i, k]) * precip_frac[i, k] +
                           nprci[i, k] * icldm[i, k] - nscng[i, k] * lcldm[i, k]
            ngtend[i, k] = ngtend[i, k] + nscng[i, k] * lcldm[i, k] +
                           (ngracs[i, k] + nnuccr[i, k]) * precip_frac[i, k]
        else
            nstend[i, k] = nstend[i, k] +
                           (nsubs[i, k] + nsagg[i, k] + nnuccr[i, k]) * precip_frac[i, k] +
                           nprci[i, k] * icldm[i, k]
        end
        nrtend[i, k] = nrtend[i, k] + nprc[i, k] * lcldm[i, k] +
                       (nsubr[i, k] - npracs[i, k] - nnuccr[i, k] - nnuccri[i, k] +
                        nragg[i, k] - npracg[i, k] - ngracs[i, k]) * precip_frac[i, k]
        qrout[i, k] = qr[i, k]
        nrout[i, k] = nr[i, k] * rho[i, k]
        qsout[i, k] = qs[i, k]
        nsout[i, k] = ns[i, k] * rho[i, k]
        qgout[i, k] = qg[i, k]
        ngout[i, k] = ng[i, k] * rho[i, k]
    end

    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qric[i, k], nric[i, k])
        lamr[i, k], n0r[i, k], nric[i, k] = rain.slope, rain.intercept, rain.number
    end
    for k in 1:levels, i in 1:columns
        rercld[i, k] = add_rain_cloud_effective_radius(rercld[i, k], lamr[i, k], n0r[i, k],
            lamc[i, k], pgam[i, k], qric[i, k], qcic[i, k], ncic[i, k])
    end
    return nothing
end
