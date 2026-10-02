# 🧰🌍 Stage 1, setup (micro_pumas_v1.F90 lines 1238–1710): working copies of the state,
#    cloud fractions, air properties, saturation and the zeroing of every output.

"""
    setup_step!(work, parameters)

🧰 Stage 1. Every container except the inputs is zeroed first; the Fortran zeroes most
fields here and writes the rest before reading them.
"""
function setup_step!(work, parameters)
    (; settings, saturation) = parameters
    (; r, rhosu) = host_short_names(parameters)
    g = parameters.host.gravity
    agraupel = dense_precipitation(parameters).a
    F = work.fortran
    for container in (work.state, work.incloud, work.rates, work.process_rates,
        work.tendencies, work.diagnostics, work.scratch)
        zero_fields!(container)
    end
    (; t, q, p, pdel, cldn, liqcldf, icecldf, qsatfac) = F
    (; qc, nc, qi, ni, qr, nr, qs, ns, qg, ng) = F
    (; cldm, lcldm, icldm, qsfm, precip_frac) = F
    (; rho, dv, mu, sc, rhof, arn, asn, agn, acn, ain, ajn) = F
    (; esl, esi, qvl, qvi, relhum) = F
    columns, levels = grid_size(work)

    # 📋 copies of the input concentrations that are changed internally
    for k in 1:levels, i in 1:columns
        qc[i, k] = F.qcn[i, k]
        nc[i, k] = F.ncn[i, k]
        qi[i, k] = F.qin[i, k]
        ni[i, k] = F.nin[i, k]
        qr[i, k] = F.qrn[i, k]
        nr[i, k] = F.nrn[i, k]
        qs[i, k] = F.qsn[i, k]
        ns[i, k] = F.nsn[i, k]
        qg[i, k] = F.qgr[i, k]
        ng[i, k] = F.ngr[i, k]
    end

    # ☁️ cloud fractions, with a minimum against division by zero
    for k in 1:levels, i in 1:columns
        if settings.uniform_subcolumns
            lcldm[i, k] = qc[i, k] >= qsmall ? 1.0 : mincld
            icldm[i, k] = qi[i, k] >= qsmall ? 1.0 : mincld
            cldm[i, k] = max(icldm[i, k], lcldm[i, k])
            qsfm[i, k] = 1.0
        else
            cldm[i, k] = max(cldn[i, k], mincld)
            lcldm[i, k] = max(liqcldf[i, k], mincld)
            icldm[i, k] = max(icecldf[i, k], mincld)
            qsfm[i, k] = qsatfac[i, k]
        end
    end

    # 🌬️ air properties and density-corrected fall speed prefactors
    for k in 1:levels, i in 1:columns
        rho[i, k] = p[i, k] / (r * t[i, k])
        dv[i, k] = 8.794e-5 * t[i, k]^1.81 / p[i, k]
        mu[i, k] = 1.496e-6 * t[i, k]^1.5 / (t[i, k] + 120.0)
        sc[i, k] = mu[i, k] / (rho[i, k] * dv[i, k])
        rhof[i, k] = (rhosu / rho[i, k])^0.54
        arn[i, k] = fall_speed_law.rain.a * rhof[i, k]
        asn[i, k] = fall_speed_law.snow.a * rhof[i, k]
        if settings.do_hail || settings.do_graupel
            agn[i, k] = agraupel * rhof[i, k]
        end
        acn[i, k] = g * rhow / (18.0 * mu[i, k])
        ain[i, k] = fall_speed_law.cloud_ice.a * (rhosu / rho[i, k])^0.35
        ajn[i, k] = fall_speed_law.small_cloud_ice.a * (rhosu / rho[i, k])^0.35
    end

    # 💧 saturation over water and ice, scaled for subgrid ice cloud below freezing
    tmelt = parameters.host.melting_point
    for k in 1:levels, i in 1:columns
        water = saturation_over_water(saturation, t[i, k], p[i, k])
        ice = saturation_over_ice(saturation, t[i, k], p[i, k])
        esl[i, k], qvl[i, k] = water.vapor_pressure, water.specific_humidity
        esi[i, k], qvi[i, k] = ice.vapor_pressure, ice.specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if t[i, k] >= tmelt
            esi[i, k] = esl[i, k]
            qvi[i, k] = qvl[i, k]
        else
            qvi[i, k] = qsfm[i, k] * qvi[i, k]
            esi[i, k] = qsfm[i, k] * esi[i, k]
            qvl[i, k] = qsfm[i, k] * qvl[i, k]
            esl[i, k] = qsfm[i, k] * esl[i, k]
        end
        relhum[i, k] = q[i, k] / max(qvl[i, k], qsmall)
    end

    # 🧹 outputs that do not start at zero
    for k in 1:levels, i in 1:columns
        precip_frac[i, k] = mincld
        F.qcrat[i, k] = 1.0
        F.effc[i, k] = 10.0
        F.effc_fn[i, k] = 10.0
        F.effi[i, k] = 25.0
        F.effi[i, k] = F.effi[i, k] * settings.ice_effective_radius_factor
        F.deffi[i, k] = 50.0
        F.refl[i, k] = -9999.0
        F.refl10cm[i, k] = -9999.0
        F.pdel_inv[i, k] = 1.0 / pdel[i, k]
    end
    return nothing
end
