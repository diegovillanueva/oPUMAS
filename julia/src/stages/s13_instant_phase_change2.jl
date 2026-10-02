# 🌡️🧊 Stage 13, instantaneous processes after sedimentation (micro_pumas_v1.F90 lines
#    3733–4020): snow and graupel melting, rain freezing, cloud ice melting, homogeneous
#    droplet freezing, the ice number limit and the removal of supersaturation.

"""
    instant_phase_change_after_sedimentation!(work, parameters)

🌡️ Stage 13. Each instantaneous change is limited so that it does not push the
end-of-step temperature across its threshold; tendencies and `tlat` are updated.
"""
function instant_phase_change_after_sedimentation!(work, parameters)
    (; settings, saturation, hydrometeors) = parameters
    (; cpp, rv, xxlv, xxls, xlf, tmelt, snowmelt, rainfrze) = host_short_names(parameters)
    xxlv_squared = xxlv^2
    F = work.fortran
    (; t, q, p, qc, nc, qi, ni, qr, nr, qs, ns, qg, ng, rho, icldm, precip_frac) = F
    (; tlat, qvlat, qctend, nctend, qitend, nitend, qrtend, nrtend, qstend, nstend) = F
    (; qgtend, ngtend, cmeout, lamr, dum_2D, ttmpA, qvnA) = F
    (; dumc, dumnc, dumi, dumni, dumr, dumnr, dums, dumns, dumg, dumng) = F
    deltat = work.inputs.time_step
    rdeltat = 1.0 / deltat
    Dcs = settings.ice_to_snow_size
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        dumc[i, k] = max(qc[i, k] + qctend[i, k] * deltat, 0.0)
        dumi[i, k] = max(qi[i, k] + qitend[i, k] * deltat, 0.0)
        dumnc[i, k] = max(nc[i, k] + nctend[i, k] * deltat, 0.0)
        dumni[i, k] = max(ni[i, k] + nitend[i, k] * deltat, 0.0)
        dumr[i, k] = max(qr[i, k] + qrtend[i, k] * deltat, 0.0)
        dumnr[i, k] = max(nr[i, k] + nrtend[i, k] * deltat, 0.0)
        dums[i, k] = max(qs[i, k] + qstend[i, k] * deltat, 0.0)
        dumns[i, k] = max(ns[i, k] + nstend[i, k] * deltat, 0.0)
        dumg[i, k] = max(qg[i, k] + qgtend[i, k] * deltat, 0.0)
        dumng[i, k] = max(ng[i, k] + ngtend[i, k] * deltat, 0.0)
        if settings.constant_droplet_number
            dumnc[i, k] = settings.droplet_number / rho[i, k] * F.lcldm[i, k]
        end
        if settings.constant_ice_number
            dumni[i, k] = settings.ice_number / rho[i, k] * icldm[i, k]
        end
        if settings.constant_graupel_number
            dumng[i, k] = settings.graupel_number / rho[i, k] * precip_frac[i, k]
        end
        if dumc[i, k] < qsmall
            dumnc[i, k] = 0.0
        end
        if dumi[i, k] < qsmall
            dumni[i, k] = 0.0
        end
        if dumr[i, k] < qsmall
            dumnr[i, k] = 0.0
        end
        if dums[i, k] < qsmall
            dumns[i, k] = 0.0
        end
        if dumg[i, k] < qsmall
            dumng[i, k] = 0.0
        end

        # ❄️💧 melting of snow at +2 °C
        if t[i, k] + tlat[i, k] / cpp * deltat > snowmelt
            if dums[i, k] > 0.0
                dum = -xlf / cpp * dums[i, k]
                if t[i, k] + tlat[i, k] / cpp * deltat + dum < snowmelt
                    dum = (t[i, k] + tlat[i, k] / cpp * deltat - snowmelt) * cpp / xlf
                    dum = dum / dums[i, k]
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else
                    dum = 1.0
                end
                qstend[i, k] = qstend[i, k] - dum * dums[i, k] * rdeltat
                nstend[i, k] = nstend[i, k] - dum * dumns[i, k] * rdeltat
                qrtend[i, k] = qrtend[i, k] + dum * dums[i, k] * rdeltat
                nrtend[i, k] = nrtend[i, k] + dum * dumns[i, k] * rdeltat
                dum1 = -xlf * dum * dums[i, k] * rdeltat
                tlat[i, k] = tlat[i, k] + dum1
                F.meltsdttot[i, k] = F.meltsdttot[i, k] + dum1
            end
        end

        # 🟤💧 melting of graupel at +2 °C
        if t[i, k] + tlat[i, k] / cpp * deltat > snowmelt
            if dumg[i, k] > 0.0
                dum = -xlf / cpp * dumg[i, k]
                if t[i, k] + tlat[i, k] / cpp * deltat + dum < snowmelt
                    dum = (t[i, k] + tlat[i, k] / cpp * deltat - snowmelt) * cpp / xlf
                    dum = dum / dumg[i, k]
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else
                    dum = 1.0
                end
                qgtend[i, k] = qgtend[i, k] - dum * dumg[i, k] * rdeltat
                ngtend[i, k] = ngtend[i, k] - dum * dumng[i, k] * rdeltat
                qrtend[i, k] = qrtend[i, k] + dum * dumg[i, k] * rdeltat
                nrtend[i, k] = nrtend[i, k] + dum * dumng[i, k] * rdeltat
                dum1 = -xlf * dum * dumg[i, k] * rdeltat
                tlat[i, k] = tlat[i, k] + dum1
                F.meltsdttot[i, k] = F.meltsdttot[i, k] + dum1
                if settings.constant_graupel_number && (ng[i, k] + ngtend[i, k] * deltat) < 0.0
                    ngtend[i, k] = -ng[i, k] * rdeltat
                end
            end
        end
    end

    # 🌧️ the rain size uses the grid-mean rain (FINDINGS.md)
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, dumr[i, k], dumnr[i, k])
        lamr[i, k], dumnr[i, k] = rain.slope, rain.number
    end

    # 🌧️🧊 freezing of rain at -5 °C
    for k in 1:levels, i in 1:columns
        if t[i, k] + tlat[i, k] / cpp * deltat < rainfrze
            if dumr[i, k] > 0.0
                dum = xlf / cpp * dumr[i, k]
                if t[i, k] + tlat[i, k] / cpp * deltat + dum > rainfrze
                    dum = -(t[i, k] + tlat[i, k] / cpp * deltat - rainfrze) * cpp / xlf
                    dum = dum / dumr[i, k]
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else
                    dum = 1.0
                end
                qrtend[i, k] = qrtend[i, k] - dum * dumr[i, k] * rdeltat
                nrtend[i, k] = nrtend[i, k] - dum * dumnr[i, k] * rdeltat
                if lamr[i, k] < 1.0 / Dcs
                    if settings.do_hail || settings.do_graupel
                        qgtend[i, k] = qgtend[i, k] + dum * dumr[i, k] * rdeltat
                        ngtend[i, k] = ngtend[i, k] + dum * dumnr[i, k] * rdeltat
                    else
                        qstend[i, k] = qstend[i, k] + dum * dumr[i, k] * rdeltat
                        nstend[i, k] = nstend[i, k] + dum * dumnr[i, k] * rdeltat
                    end
                else
                    qitend[i, k] = qitend[i, k] + dum * dumr[i, k] * rdeltat
                    nitend[i, k] = nitend[i, k] + dum * dumnr[i, k] * rdeltat
                end
                dum1 = xlf * dum * dumr[i, k] * rdeltat
                F.frzrdttot[i, k] = F.frzrdttot[i, k] + dum1
                tlat[i, k] = tlat[i, k] + dum1
            end
        end
    end

    settings.prognostic_cloud_ice || return nothing
    homog_size = settings.homogeneous_freezing_size
    max_nicons = settings.maximum_ice_number
    for k in 1:levels, i in 1:columns
        # 🧊💧 melting of cloud ice
        if t[i, k] + tlat[i, k] / cpp * deltat > tmelt
            if dumi[i, k] > 0.0
                dum = -dumi[i, k] * xlf / cpp
                if t[i, k] + tlat[i, k] / cpp * deltat + dum < tmelt
                    dum = (t[i, k] + tlat[i, k] / cpp * deltat - tmelt) * cpp / xlf
                    dum = dum / dumi[i, k]
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else
                    dum = 1.0
                end
                qctend[i, k] = qctend[i, k] + dum * dumi[i, k] * rdeltat
                F.melttot[i, k] = dum * dumi[i, k] * rdeltat
                F.nmelttot[i, k] = 3.0 * dum * dumi[i, k] * rdeltat /
                                   (4.0 * π * 5.12e-16 * rhow)
                nctend[i, k] = nctend[i, k] + F.nmelttot[i, k]
                qitend[i, k] = ((1.0 - dum) * dumi[i, k] - qi[i, k]) * rdeltat
                nitend[i, k] = ((1.0 - dum) * dumni[i, k] - ni[i, k]) * rdeltat
                tlat[i, k] = tlat[i, k] - xlf * dum * dumi[i, k] * rdeltat
            end
        end

        # 💧➡️🧊 homogeneous freezing of droplets at -40 °C
        if t[i, k] + tlat[i, k] / cpp * deltat < 233.15
            if dumc[i, k] > 0.0
                dum = dumc[i, k] * xlf / cpp
                if t[i, k] + tlat[i, k] / cpp * deltat + dum > 233.15
                    dum = -(t[i, k] + tlat[i, k] / cpp * deltat - 233.15) * cpp / xlf
                    dum = dum / dumc[i, k]
                    dum = max(0.0, dum)
                    dum = min(1.0, dum)
                else
                    dum = 1.0
                end
                qitend[i, k] = qitend[i, k] + dum * dumc[i, k] * rdeltat
                F.homotot[i, k] = dum * dumc[i, k] * rdeltat
                F.nhomotot[i, k] = dum * 3.0 * dumc[i, k] /
                                   (4.0 * 3.14 * homog_size^3.0 * 500.0) * rdeltat
                nitend[i, k] = nitend[i, k] + F.nhomotot[i, k]
                qctend[i, k] = ((1.0 - dum) * dumc[i, k] - qc[i, k]) * rdeltat
                nctend[i, k] = ((1.0 - dum) * dumnc[i, k] - nc[i, k]) * rdeltat
                tlat[i, k] = tlat[i, k] + xlf * dum * dumc[i, k] * rdeltat
            end
        end

        # 🔢 ice number limiter
        if nitend[i, k] > 0.0 &&
           ni[i, k] + nitend[i, k] * deltat > max_nicons * icldm[i, k] / rho[i, k]
            nitend[i, k] = max(0.0, (max_nicons * icldm[i, k] / rho[i, k] - ni[i, k]) / deltat)
        end

        dum_2D[i, k] = q[i, k] + qvlat[i, k] * deltat
        ttmpA[i, k] = t[i, k] + tlat[i, k] / cpp * deltat
    end

    for k in 1:levels, i in 1:columns
        qvnA[i, k] = saturation_over_water(saturation, ttmpA[i, k], p[i, k]).specific_humidity
    end

    # 🌫️ remove the supersaturation left by the sum of all processes
    for k in 1:levels, i in 1:columns
        if dum_2D[i, k] > qvnA[i, k] && qvnA[i, k] > 0 && settings.remove_supersaturation
            dum = (dum_2D[i, k] - qvnA[i, k]) /
                  (1.0 + xxlv_squared * qvnA[i, k] / (cpp * rv * (ttmpA[i, k] * ttmpA[i, k]))) *
                  rdeltat
            cmeout[i, k] = cmeout[i, k] + dum
            if ttmpA[i, k] > 268.15
                dum1 = 0.0
            elseif ttmpA[i, k] < 238.15
                dum1 = 1.0
            else
                dum1 = (268.15 - ttmpA[i, k]) / 30.0
            end
            latent = xxls * dum1 + xxlv * (1.0 - dum1)
            dum = (dum_2D[i, k] - qvnA[i, k]) /
                  (1.0 + latent * latent * qvnA[i, k] / (cpp * rv * (ttmpA[i, k] * ttmpA[i, k]))) *
                  rdeltat
            qctend[i, k] = qctend[i, k] + dum * (1.0 - dum1)
            F.qcrestot[i, k] = dum * (1.0 - dum1)
            qitend[i, k] = qitend[i, k] + dum * dum1
            F.qirestot[i, k] = dum * dum1
            qvlat[i, k] = qvlat[i, k] - dum
            F.qvres[i, k] = -dum
            tlat[i, k] = tlat[i, k] + dum * (1.0 - dum1) * xxlv + dum * dum1 * xxls
        end
    end
    return nothing
end
