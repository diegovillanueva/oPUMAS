# 🪂📏 Stage 11, fall speeds for sedimentation (micro_pumas_v1.F90 lines 3222–3505): the
#    end-of-step in-cloud values give new size distributions, then the mass and number
#    weighted fall speeds as pressure fluxes g ρ V. The dummies become grid means.

"""
    sedimentation_fall_speeds!(work, parameters)

🪂 Stage 11. Fills `dumc`…`dumng` (end-of-step in-cloud, then grid-mean values),
`fc`…`fng` (g ρ V) and the fall speed process rates, and `zint` for the implicit
scheme.
"""
function sedimentation_fall_speeds!(work, parameters)
    (; settings, hydrometeors, gamma) = parameters
    (; g, r) = host_short_names(parameters)
    dense = dense_precipitation(parameters)
    F = work.fortran
    (; t, pint, qc, nc, qi, ni, qr, nr, qs, ns, qg, ng, rho, rhof) = F
    (; acn, ain, ajn, arn, asn, agn, lcldm, icldm, precip_frac) = F
    (; qctend, nctend, qitend, nitend, qrtend, nrtend, qstend, nstend, qgtend, ngtend) = F
    (; dumc, dumnc, dumi, dumni, dumr, dumnr, dums, dumns, dumg, dumng) = F
    (; fc, fnc, fi, fni, fr, fnr, fs, fns, fg, fng) = F
    (; lami, pgam, lamc, lamr, lams, lamg, unr, uns, ung, zint) = F
    (; vtrmc, vtrmi, umr, ums, umg) = F
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        dumc[i, k] = (qc[i, k] + qctend[i, k] * deltat) / lcldm[i, k]
        dumi[i, k] = (qi[i, k] + qitend[i, k] * deltat) / icldm[i, k]
        dumnc[i, k] = max((nc[i, k] + nctend[i, k] * deltat) / lcldm[i, k], 0.0)
        dumni[i, k] = max((ni[i, k] + nitend[i, k] * deltat) / icldm[i, k], 0.0)
        dumr[i, k] = (qr[i, k] + qrtend[i, k] * deltat) / precip_frac[i, k]
        dumnr[i, k] = max((nr[i, k] + nrtend[i, k] * deltat) / precip_frac[i, k], 0.0)
        dums[i, k] = (qs[i, k] + qstend[i, k] * deltat) / precip_frac[i, k]
        dumns[i, k] = max((ns[i, k] + nstend[i, k] * deltat) / precip_frac[i, k], 0.0)
        dumg[i, k] = (qg[i, k] + qgtend[i, k] * deltat) / precip_frac[i, k]
        dumng[i, k] = max((ng[i, k] + ngtend[i, k] * deltat) / precip_frac[i, k], 0.0)
        # 🔢 constant numbers; the constant rain and snow numbers (nrcons, nscons) are
        #    not settings of oPUMAS.jl yet, so their branches are left out.
        if settings.constant_graupel_number
            dumng[i, k] = settings.graupel_number / rho[i, k]
        end
        if settings.constant_droplet_number
            dumnc[i, k] = settings.droplet_number / rho[i, k]
        end
        if settings.constant_ice_number
            dumni[i, k] = settings.ice_number / rho[i, k]
        end
    end

    for k in 1:levels, i in 1:columns
        ice = size_distribution(hydrometeors.cloud_ice, dumi[i, k], dumni[i, k])
        lami[i, k], dumni[i, k] = ice.slope, ice.number
    end
    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, dumc[i, k],
            dumnc[i, k], rho[i, k])
        pgam[i, k], lamc[i, k], dumnc[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, dumr[i, k], dumnr[i, k])
        lamr[i, k], dumnr[i, k] = rain.slope, rain.number
    end
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, dums[i, k], dumns[i, k])
        lams[i, k], dumns[i, k] = snow.slope, snow.number
    end
    if settings.do_hail || settings.do_graupel
        dense_props = settings.do_hail ? hydrometeors.hail : hydrometeors.graupel
        for k in 1:levels, i in 1:columns
            graupel = size_distribution(dense_props, dumg[i, k], dumng[i, k])
            lamg[i, k], dumng[i, k] = graupel.slope, graupel.number
        end
    end

    if settings.implicit_sedimentation
        for i in 1:columns
            zint[i, levels+1] = 0.0
            for k in levels:-1:1
                H = r * t[i, k] / g * log(pint[i, k+1] / pint[i, k])
                zint[i, k] = zint[i, k+1] + H
            end
        end
    end

    # 💧 droplets
    bc = fall_speed_law.cloud_liquid.b
    for k in 1:levels, i in 1:columns
        if dumc[i, k] >= qsmall
            dum1 = 4.0 + bc + pgam[i, k]
            dum2 = pgam[i, k] + 4.0
            vtrmc[i, k] = acn[i, k] * gamma_function(dum1) /
                          (lamc[i, k]^bc * gamma_function(dum2))
            if settings.ifs_snow_fall_speed
                fc[i, k] = 0.0
                fnc[i, k] = 0.0
            else
                dum3 = 1.0 + bc + pgam[i, k]
                dum4 = pgam[i, k] + 1.0
                fc[i, k] = g * rho[i, k] * vtrmc[i, k]
                fnc[i, k] = g * rho[i, k] * acn[i, k] * gamma_function(dum3) /
                            (lamc[i, k]^bc * gamma_function(dum4))
            end
        else
            fc[i, k] = 0.0
            fnc[i, k] = 0.0
        end
        dumc[i, k] = (qc[i, k] + qctend[i, k] * deltat)
        dumnc[i, k] = max((nc[i, k] + nctend[i, k] * deltat), 0.0)
        if dumc[i, k] < qsmall
            dumnc[i, k] = 0.0
        end
    end

    # 🧊 cloud ice, blended with small-ice fall speeds below 20 μm
    bi, bj = fall_speed_law.cloud_ice.b, fall_speed_law.small_cloud_ice.b
    vtrmi_factor = settings.ice_fall_speed_factor
    for k in 1:levels, i in 1:columns
        if dumi[i, k] >= qsmall
            vtrmi[i, k] = min(ain[i, k] * gamma.ice_plus4 / (6.0 * lami[i, k]^bi),
                1.2 * rhof[i, k])
            vtrmi[i, k] = vtrmi[i, k] * vtrmi_factor
            fi[i, k] = g * rho[i, k] * vtrmi[i, k]
            fni[i, k] = g * rho[i, k] *
                        min(ain[i, k] * gamma.ice_plus1 / lami[i, k]^bi, 1.2 * rhof[i, k])
            irad = 1.5 / lami[i, k] * 1e6
            ifrac = min(1.0, max(0.0, (irad - 18.0) / 2.0))
            if ifrac < 1.0
                # 🔁 the factor is applied again to the blend (FINDINGS.md)
                vtrmi[i, k] = ifrac * vtrmi[i, k] +
                              (1.0 - ifrac) *
                              min(ajn[i, k] * gamma.small_ice_plus4 / (6.0 * lami[i, k]^bj),
                    1.2 * rhof[i, k])
                vtrmi[i, k] = vtrmi[i, k] * vtrmi_factor
                fi[i, k] = g * rho[i, k] * vtrmi[i, k]
                fni[i, k] = ifrac * fni[i, k] +
                            (1.0 - ifrac) * g * rho[i, k] *
                            min(ajn[i, k] * gamma.small_ice_plus1 / lami[i, k]^bj,
                    1.2 * rhof[i, k])
            end
            if settings.ifs_snow_fall_speed
                fi[i, k] = g * rho[i, k] * 0.1
                fni[i, k] = g * rho[i, k] * 0.1
            end
        else
            fi[i, k] = 0.0
            fni[i, k] = 0.0
        end
        dumi[i, k] = (qi[i, k] + qitend[i, k] * deltat)
        dumni[i, k] = max((ni[i, k] + nitend[i, k] * deltat), 0.0)
        if dumi[i, k] < qsmall
            dumni[i, k] = 0.0
        end
    end

    # 🌧️ rain; the fall correction starts at k = 3 as in the Fortran (FINDINGS.md)
    br = fall_speed_law.rain.b
    for i in 1:columns, k in 1:levels
        if lamr[i, k] >= qsmall
            qtmp = lamr[i, k]^br
            unr[i, k] = min(arn[i, k] * gamma.rain_plus1 / qtmp, 9.1 * rhof[i, k])
            fnr[i, k] = g * rho[i, k] * unr[i, k]
            umr[i, k] = min(arn[i, k] * gamma.rain_plus4 / (6.0 * qtmp), 9.1 * rhof[i, k])
            fr[i, k] = g * rho[i, k] * umr[i, k]
        else
            fr[i, k] = 0.0
            fnr[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fr[i, k] < 1.e-10
            fr[i, k] = fr[i, k-1]
            fnr[i, k] = fnr[i, k-1]
        end
        dumr[i, k] = (qr[i, k] + qrtend[i, k] * deltat)
        dumnr[i, k] = max((nr[i, k] + nrtend[i, k] * deltat), 0.0)
        if dumr[i, k] < qsmall
            dumnr[i, k] = 0.0
        end
    end

    # ❄️ snow
    bs = fall_speed_law.snow.b
    for i in 1:columns, k in 1:levels
        if lams[i, k] >= qsmall
            qtmp = lams[i, k]^bs
            ums[i, k] = min(asn[i, k] * gamma.snow_plus4 / (6.0 * qtmp), 1.2 * rhof[i, k])
            ums[i, k] = ums[i, k] * settings.snow_fall_speed_factor
            fs[i, k] = g * rho[i, k] * ums[i, k]
            uns[i, k] = min(asn[i, k] * gamma.snow_plus1 / qtmp, 1.2 * rhof[i, k])
            fns[i, k] = g * rho[i, k] * uns[i, k]
            if settings.ifs_snow_fall_speed
                ums[i, k] = 1.0
                uns[i, k] = 1.0
            end
        else
            fs[i, k] = 0.0
            fns[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fs[i, k] < 1.e-10
            fs[i, k] = fs[i, k-1]
            fns[i, k] = fns[i, k-1]
        end
        dums[i, k] = (qs[i, k] + qstend[i, k] * deltat)
        dumns[i, k] = max((ns[i, k] + nstend[i, k] * deltat), 0.0)
        if dums[i, k] < qsmall
            dumns[i, k] = 0.0
        end
    end

    # 🟤 graupel or hail
    for i in 1:columns, k in 1:levels
        if lamg[i, k] >= qsmall
            qtmp = lamg[i, k]^dense.b
            umg[i, k] = min(agn[i, k] * gamma.graupel_plus4 / (6.0 * qtmp),
                20.0 * rhof[i, k])
            fg[i, k] = g * rho[i, k] * umg[i, k]
            ung[i, k] = min(agn[i, k] * gamma.graupel_plus1 / qtmp, 20.0 * rhof[i, k])
            fng[i, k] = g * rho[i, k] * ung[i, k]
        else
            fg[i, k] = 0.0
            fng[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fg[i, k] < 1.e-10
            fg[i, k] = fg[i, k-1]
            fng[i, k] = fng[i, k-1]
        end
        dumg[i, k] = (qg[i, k] + qgtend[i, k] * deltat)
        dumng[i, k] = max((ng[i, k] + ngtend[i, k] * deltat), 0.0)
        if dumg[i, k] < qsmall
            dumng[i, k] = 0.0
        end
    end
    return nothing
end
