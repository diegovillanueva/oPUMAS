# 🧰🌍 Stage 1, setup (micro_pumas_v1.F90 lines 1238–1710): working copies of the state,
#    cloud fractions, air properties, saturation and the zeroing of every output.

"""
    setup_step!(work, parameters)

🧰 Stage 1. Every container except the inputs is zeroed first; the Fortran zeroes most
fields here and writes the rest before reading them.
"""
function setup_step!(work, parameters)
    (; settings, saturation) = parameters
    (; Rd, ρ_ref) = host_constants(parameters)
    g = parameters.host.gravity
    agraupel = dense_precipitation(parameters).a
    (; inputs, state, incloud, diagnostics, scratch) = work
    qc_in = inputs.cloud_liquid.q  # qcn
    nc_in = inputs.cloud_liquid.N  # ncn
    qi_in = inputs.cloud_ice.q  # qin
    ni_in = inputs.cloud_ice.N  # nin
    qr_in = inputs.rain.q  # qrn
    nr_in = inputs.rain.N  # nrn
    qs_in = inputs.snow.q  # qsn
    ns_in = inputs.snow.N  # nsn
    qg_in = inputs.graupel.q  # qgr
    ng_in = inputs.graupel.N  # ngr
    qc_limiter_ratio = diagnostics.cloud_liquid.limiter_ratio  # qcrat
    radius_cloud_liquid = diagnostics.effective_radius.cloud_liquid  # effc
    radius_cloud_liquid_fixed = diagnostics.effective_radius.cloud_liquid_fixed_number  # effc_fn
    radius_cloud_ice = diagnostics.effective_radius.cloud_ice  # effi
    diameter_cloud_ice = diagnostics.effective_diameter.cloud_ice  # deffi
    dbz = diagnostics.radar.reflectivity_dbz  # refl
    dbz_10cm = diagnostics.radar.reflectivity_10cm_dbz  # refl10cm
    inv_Δp = scratch.inverse_pressure_thickness  # pdel_inv
    for container in (work.state, work.incloud, work.rates, work.process_rates,
        work.tendencies, work.diagnostics, work.scratch)
        zero_fields!(container)
    end
    T = inputs.T  # t
    qv = inputs.q_vapor  # q
    p = inputs.pressure
    Δp = inputs.pressure_thickness  # pdel
    cloud_fraction_in = inputs.cloud_fraction  # cldn
    liquid_fraction_in = inputs.liquid_cloud_fraction  # liqcldf
    ice_fraction_in = inputs.ice_cloud_fraction  # icecldf
    qsat_factor_in = inputs.saturation_scaling  # qsatfac
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    qi = state.cloud_ice.q
    ni = state.cloud_ice.N
    qr = state.rain.q
    nr = state.rain.N
    qs = state.snow.q
    ns = state.snow.N
    qg = state.graupel.q
    ng = state.graupel.N
    cloud_fraction = incloud.cloud_fraction  # cldm
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    qsat_factor = incloud.saturation_scaling  # qsfm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    ρ = state.air.density  # rho
    D_v = state.air.diffusivity  # dv
    η = state.air.viscosity  # mu
    Sc = state.air.schmidt  # sc
    ρ_correction = state.air.density_correction  # rhof
    a_rain = state.air.fall_speed_prefactor.rain  # arn
    a_snow = state.air.fall_speed_prefactor.snow  # asn
    a_graupel = state.air.fall_speed_prefactor.graupel  # agn
    a_cloud_liquid = state.air.fall_speed_prefactor.cloud_liquid  # acn
    a_cloud_ice = state.air.fall_speed_prefactor.cloud_ice  # ain
    a_small_ice = state.air.fall_speed_prefactor.small_cloud_ice  # ajn
    e_sat_liquid = state.air.vapor_pressure_liquid  # esl
    e_sat_ice = state.air.vapor_pressure_ice  # esi
    qsat_liquid = state.air.q_sat_liquid  # qvl
    qsat_ice = state.air.q_sat_ice  # qvi
    rh = state.air.relative_humidity  # relhum
    columns, levels = grid_size(work)

    # 📋 copies of the input concentrations that are changed internally
    for k in 1:levels, i in 1:columns
        qc[i, k] = qc_in[i, k]
        nc[i, k] = nc_in[i, k]
        qi[i, k] = qi_in[i, k]
        ni[i, k] = ni_in[i, k]
        qr[i, k] = qr_in[i, k]
        nr[i, k] = nr_in[i, k]
        qs[i, k] = qs_in[i, k]
        ns[i, k] = ns_in[i, k]
        qg[i, k] = qg_in[i, k]
        ng[i, k] = ng_in[i, k]
    end

    # ☁️ cloud fractions, with a minimum against division by zero
    for k in 1:levels, i in 1:columns
        if settings.uniform_subcolumns
            liquid_fraction[i, k] = qc[i, k] >= smallest_mixing_ratio ? 1.0 : minimum_cloud_fraction
            ice_fraction[i, k] = qi[i, k] >= smallest_mixing_ratio ? 1.0 : minimum_cloud_fraction
            cloud_fraction[i, k] = max(ice_fraction[i, k], liquid_fraction[i, k])
            qsat_factor[i, k] = 1.0
        else
            cloud_fraction[i, k] = max(cloud_fraction_in[i, k], minimum_cloud_fraction)
            liquid_fraction[i, k] = max(liquid_fraction_in[i, k], minimum_cloud_fraction)
            ice_fraction[i, k] = max(ice_fraction_in[i, k], minimum_cloud_fraction)
            qsat_factor[i, k] = qsat_factor_in[i, k]
        end
    end

    # 🌬️ air properties and density-corrected fall speed prefactors
    for k in 1:levels, i in 1:columns
        ρ[i, k] = p[i, k] / (Rd * T[i, k])
        D_v[i, k] = 8.794e-5 * T[i, k]^1.81 / p[i, k]
        η[i, k] = 1.496e-6 * T[i, k]^1.5 / (T[i, k] + 120.0)
        Sc[i, k] = η[i, k] / (ρ[i, k] * D_v[i, k])
        ρ_correction[i, k] = (ρ_ref / ρ[i, k])^0.54
        a_rain[i, k] = fall_speed_law.rain.a * ρ_correction[i, k]
        a_snow[i, k] = fall_speed_law.snow.a * ρ_correction[i, k]
        if settings.do_hail || settings.do_graupel
            a_graupel[i, k] = agraupel * ρ_correction[i, k]
        end
        a_cloud_liquid[i, k] = g * ρ_water / (18.0 * η[i, k])
        a_cloud_ice[i, k] = fall_speed_law.cloud_ice.a * (ρ_ref / ρ[i, k])^0.35
        a_small_ice[i, k] = fall_speed_law.small_cloud_ice.a * (ρ_ref / ρ[i, k])^0.35
    end

    # 💧 saturation over water and ice, scaled for subgrid ice cloud below freezing
    T₀ = parameters.host.melting_point
    for k in 1:levels, i in 1:columns
        water = saturation_over_water(saturation, T[i, k], p[i, k])
        ice = saturation_over_ice(saturation, T[i, k], p[i, k])
        e_sat_liquid[i, k], qsat_liquid[i, k] = water.vapor_pressure, water.specific_humidity
        e_sat_ice[i, k], qsat_ice[i, k] = ice.vapor_pressure, ice.specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if T[i, k] >= T₀
            e_sat_ice[i, k] = e_sat_liquid[i, k]
            qsat_ice[i, k] = qsat_liquid[i, k]
        else
            qsat_ice[i, k] = qsat_factor[i, k] * qsat_ice[i, k]
            e_sat_ice[i, k] = qsat_factor[i, k] * e_sat_ice[i, k]
            qsat_liquid[i, k] = qsat_factor[i, k] * qsat_liquid[i, k]
            e_sat_liquid[i, k] = qsat_factor[i, k] * e_sat_liquid[i, k]
        end
        rh[i, k] = qv[i, k] / max(qsat_liquid[i, k], smallest_mixing_ratio)
    end

    # 🧹 outputs that do not start at zero
    for k in 1:levels, i in 1:columns
        precip_fraction[i, k] = minimum_cloud_fraction
        qc_limiter_ratio[i, k] = 1.0
        radius_cloud_liquid[i, k] = 10.0
        radius_cloud_liquid_fixed[i, k] = 10.0
        radius_cloud_ice[i, k] = 25.0
        radius_cloud_ice[i, k] = radius_cloud_ice[i, k] * settings.ice_effective_radius_factor
        diameter_cloud_ice[i, k] = 50.0
        dbz[i, k] = -9999.0
        dbz_10cm[i, k] = -9999.0
        inv_Δp[i, k] = 1.0 / Δp[i, k]
    end
    return nothing
end
