# 📏🌈 Stage 14, effective radii for radiation and final number fixes
#    (micro_pumas_v1.F90 lines 4022–4307): end-of-step in-cloud values give the size
#    distributions, the effective radii, and number tendencies that keep sizes in range.

"""
    effective_radius!(work, parameters)

📏 Stage 14. Fills `effc`, `effc_fn`, `effi`, `deffi`, `sadice`, `sadsnow`, `lamcrad`,
`pgamrad`, adjusts the number tendencies, and adds the tendencies to `qc` and `qi`
(they are only used for output from here on).
"""
function effective_radius!(work, parameters)
    (; settings, hydrometeors) = parameters
    F = work.fortran
    (; qc, nc, qi, ni, qr, nr, qs, ns, qg, ng, rho, re_ice, lcldm, icldm, precip_frac) = F
    (; qctend, nctend, qitend, nitend, qrtend, nrtend, qstend, nstend, qgtend, ngtend) = F
    (; dumc, dumnc, dumi, dumni, dumr, dumnr, dums, dumns, dumg, dumng, dum_2D) = F
    (; lami, pgam, lamc, lamr, n0r, lams, lamg, dumni0A2D, dumns0A2D) = F
    (; effc, effc_fn, effi, deffi, sadice, sadsnow, lamcrad, pgamrad, nrout) = F
    deltat = work.inputs.time_step
    rdeltat = 1.0 / deltat
    effi_factor = settings.ice_effective_radius_factor
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        dumc[i, k] = max(qc[i, k] + qctend[i, k] * deltat, 0.0) / lcldm[i, k]
        dumi[i, k] = max(qi[i, k] + qitend[i, k] * deltat, 0.0) / icldm[i, k]
        dumnc[i, k] = max(nc[i, k] + nctend[i, k] * deltat, 0.0) / lcldm[i, k]
        dumni[i, k] = max(ni[i, k] + nitend[i, k] * deltat, 0.0) / icldm[i, k]
        dumr[i, k] = max(qr[i, k] + qrtend[i, k] * deltat, 0.0) / precip_frac[i, k]
        dumnr[i, k] = max(nr[i, k] + nrtend[i, k] * deltat, 0.0) / precip_frac[i, k]
        dums[i, k] = max(qs[i, k] + qstend[i, k] * deltat, 0.0) / precip_frac[i, k]
        dumns[i, k] = max(ns[i, k] + nstend[i, k] * deltat, 0.0) / precip_frac[i, k]
        # 🟤 graupel stays a grid mean here (FINDINGS.md)
        dumg[i, k] = max(qg[i, k] + qgtend[i, k] * deltat, 0.0)
        dumng[i, k] = max(ng[i, k] + ngtend[i, k] * deltat, 0.0)
        if settings.constant_droplet_number
            dumnc[i, k] = settings.droplet_number / rho[i, k]
        end
        if settings.constant_ice_number
            dumni[i, k] = settings.ice_number / rho[i, k]
        end
        if settings.constant_graupel_number
            dumng[i, k] = settings.graupel_number / rho[i, k] * precip_frac[i, k]
        end
        dumc[i, k] = min(dumc[i, k], 5.0e-3)
        dumi[i, k] = min(dumi[i, k], 5.0e-3)
        dumr[i, k] = min(dumr[i, k], 10.0e-3)
        dums[i, k] = min(dums[i, k], 10.0e-3)
        dumg[i, k] = min(dumg[i, k], 10.0e-3)
    end

    # 🧊 cloud ice
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            dum_2D[i, k] = dumni[i, k]
            ice = size_distribution(hydrometeors.cloud_ice, dumi[i, k], dumni[i, k])
            lami[i, k], dumni[i, k], dumni0A2D[i, k] = ice.slope, ice.number, ice.intercept
        end
        for k in 1:levels, i in 1:columns
            if dumi[i, k] >= qsmall
                if dumni[i, k] != dum_2D[i, k]
                    nitend[i, k] = (dumni[i, k] * icldm[i, k] - ni[i, k]) * rdeltat
                end
                effi[i, k] = 1.5 / lami[i, k] * 1.0e6
                effi[i, k] = effi[i, k] * effi_factor
                x = lami[i, k]
                sadice[i, k] = 2.0 * π * (1.0 / (x * x * x)) * dumni0A2D[i, k] * rho[i, k] *
                               1.0e-2
            else
                effi[i, k] = 25.0
                effi[i, k] = effi[i, k] * effi_factor
                sadice[i, k] = 0.0
            end
            deffi[i, k] = effi[i, k] * rhoi / rhows * 2.0
        end
    else
        for k in 1:levels, i in 1:columns
            effi[i, k] = re_ice[i, k] * 1.0e6
            effi[i, k] = effi[i, k] * effi_factor
            deffi[i, k] = effi[i, k] * 2.0
            sadice[i, k] = 4.0 * π * (effi[i, k] * effi[i, k]) * ni[i, k] * rho[i, k] * 1e-2
        end
    end

    # 💧 droplets, then again with a fixed number of 10⁸ kg⁻¹
    for k in 1:levels, i in 1:columns
        dum_2D[i, k] = dumnc[i, k]
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, dumc[i, k],
            dumnc[i, k], rho[i, k])
        pgam[i, k], lamc[i, k], dumnc[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        if dumc[i, k] >= qsmall
            if settings.constant_droplet_number
                nctend[i, k] = (settings.droplet_number / rho[i, k] * lcldm[i, k] - nc[i, k]) *
                               rdeltat
            end
            if dum_2D[i, k] != dumnc[i, k]
                nctend[i, k] = (dumnc[i, k] * lcldm[i, k] - nc[i, k]) * rdeltat
            end
            effc[i, k] = (pgam[i, k] + 3.0) / lamc[i, k] / 2.0 * 1.0e6
            lamcrad[i, k] = lamc[i, k]
            pgamrad[i, k] = pgam[i, k]
            dumnc[i, k] = 1.0e8
        end
    end
    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, dumc[i, k],
            dumnc[i, k], rho[i, k])
        pgam[i, k], lamc[i, k], dumnc[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        if dumc[i, k] >= qsmall
            effc_fn[i, k] = (pgam[i, k] + 3.0) / lamc[i, k] / 2.0 * 1.0e6
        else
            effc[i, k] = 10.0
            lamcrad[i, k] = 0.0
            pgamrad[i, k] = 0.0
            effc_fn[i, k] = 10.0
        end
        dum_2D[i, k] = dumnr[i, k]
    end

    # 🌧️ rain
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, dumr[i, k], dumnr[i, k])
        lamr[i, k], dumnr[i, k], n0r[i, k] = rain.slope, rain.number, rain.intercept
    end
    for k in 1:levels, i in 1:columns
        if dumr[i, k] >= qsmall && dum_2D[i, k] != dumnr[i, k]
            nrtend[i, k] = (dumnr[i, k] * precip_frac[i, k] - nr[i, k]) * rdeltat
        end
        dum_2D[i, k] = dumns[i, k]
    end

    # ❄️ snow
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, dums[i, k], dumns[i, k])
        lams[i, k], dumns[i, k], dumns0A2D[i, k] = snow.slope, snow.number, snow.intercept
    end
    for k in 1:levels, i in 1:columns
        if dums[i, k] >= qsmall
            if dum_2D[i, k] != dumns[i, k]
                nstend[i, k] = (dumns[i, k] * precip_frac[i, k] - ns[i, k]) * rdeltat
            end
            x = lams[i, k]
            sadsnow[i, k] = 2.0 * π * (1.0 / (x * x * x)) * dumns0A2D[i, k] * rho[i, k] *
                            1.0e-2
        end
        dum_2D[i, k] = dumng[i, k]
    end

    # 🟤 graupel or hail (hail first, then graupel, as the Fortran calls them)
    for (on, props) in ((settings.do_hail, hydrometeors.hail),
        (settings.do_graupel, hydrometeors.graupel))
        on || continue
        for k in 1:levels, i in 1:columns
            dense = size_distribution(props, dumg[i, k], dumng[i, k])
            lamg[i, k], dumng[i, k] = dense.slope, dense.number
        end
    end

    for k in 1:levels, i in 1:columns
        if dumg[i, k] >= qsmall && dum_2D[i, k] != dumng[i, k]
            F.ngtend[i, k] = (dumng[i, k] * precip_frac[i, k] - ng[i, k]) * rdeltat
        end
        # 🚫 no number without mass
        if qc[i, k] + qctend[i, k] * deltat < qsmall
            nctend[i, k] = -nc[i, k] * rdeltat
        end
        if settings.prognostic_cloud_ice && qi[i, k] + qitend[i, k] * deltat < qsmall
            nitend[i, k] = -ni[i, k] * rdeltat
        end
        if qr[i, k] + qrtend[i, k] * deltat < qsmall
            nrtend[i, k] = -nr[i, k] * rdeltat
        end
        if qs[i, k] + qstend[i, k] * deltat < qsmall
            nstend[i, k] = -ns[i, k] * rdeltat
        end
        if qg[i, k] + qgtend[i, k] * deltat < qsmall
            F.ngtend[i, k] = -ng[i, k] * rdeltat
        end
        # 📤 qc and qi are only used for output from here on
        qc[i, k] = qc[i, k] + qctend[i, k] * deltat
        qi[i, k] = qi[i, k] + qitend[i, k] * deltat
        if nrout[i, k] == 0.0
            nrout[i, k] = 1.0e-34
        end
    end
    return nothing
end
