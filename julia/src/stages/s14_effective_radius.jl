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
    (; inputs, state, incloud, tendencies, diagnostics, scratch) = work
    ng_tendency = tendencies.graupel.N  # ngtend
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
    ρ = state.air.density  # rho
    ice_radius_in = inputs.ice_effective_radius  # re_ice
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    qc_tendency = tendencies.cloud_liquid.q  # qctend
    nc_tendency = tendencies.cloud_liquid.N  # nctend
    qi_tendency = tendencies.cloud_ice.q  # qitend
    ni_tendency = tendencies.cloud_ice.N  # nitend
    qr_tendency = tendencies.rain.q  # qrtend
    nr_tendency = tendencies.rain.N  # nrtend
    qs_tendency = tendencies.snow.q  # qstend
    ns_tendency = tendencies.snow.N  # nstend
    qg_tendency = tendencies.graupel.q  # qgtend
    qc_sed = scratch.sedimenting.cloud_liquid.q  # dumc
    nc_sed = scratch.sedimenting.cloud_liquid.N  # dumnc
    qi_sed = scratch.sedimenting.cloud_ice.q  # dumi
    ni_sed = scratch.sedimenting.cloud_ice.N  # dumni
    qr_sed = scratch.sedimenting.rain.q  # dumr
    nr_sed = scratch.sedimenting.rain.N  # dumnr
    qs_sed = scratch.sedimenting.snow.q  # dums
    ns_sed = scratch.sedimenting.snow.N  # dumns
    qg_sed = scratch.sedimenting.graupel.q  # dumg
    ng_sed = scratch.sedimenting.graupel.N  # dumng
    number_before = scratch.number_before  # dum_2D
    λi = incloud.cloud_ice.λ  # lami
    μc = incloud.cloud_liquid.μ  # pgam
    λc = incloud.cloud_liquid.λ  # lamc
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    λs = incloud.snow.λ  # lams
    λg = incloud.graupel.λ  # lamg
    ni_intercept = scratch.ice_intercept  # dumni0A2D
    ns_intercept = scratch.snow_intercept  # dumns0A2D
    radius_cloud_liquid = diagnostics.effective_radius.cloud_liquid  # effc
    radius_cloud_liquid_fixed = diagnostics.effective_radius.cloud_liquid_fixed_number  # effc_fn
    radius_cloud_ice = diagnostics.effective_radius.cloud_ice  # effi
    diameter_cloud_ice = diagnostics.effective_diameter.cloud_ice  # deffi
    area_density_ice = diagnostics.surface_area_density.cloud_ice  # sadice
    area_density_snow = diagnostics.surface_area_density.snow  # sadsnow
    λc_radiation = diagnostics.radiation.droplet_slope  # lamcrad
    μc_radiation = diagnostics.radiation.droplet_shape  # pgamrad
    nr_per_volume = diagnostics.rain.N_per_volume  # nrout
    deltat = work.inputs.time_step
    rdeltat = 1.0 / deltat
    effi_factor = settings.ice_effective_radius_factor
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        qc_sed[i, k] = max(qc[i, k] + qc_tendency[i, k] * deltat, 0.0) / liquid_fraction[i, k]
        qi_sed[i, k] = max(qi[i, k] + qi_tendency[i, k] * deltat, 0.0) / ice_fraction[i, k]
        nc_sed[i, k] = max(nc[i, k] + nc_tendency[i, k] * deltat, 0.0) / liquid_fraction[i, k]
        ni_sed[i, k] = max(ni[i, k] + ni_tendency[i, k] * deltat, 0.0) / ice_fraction[i, k]
        qr_sed[i, k] = max(qr[i, k] + qr_tendency[i, k] * deltat, 0.0) / precip_fraction[i, k]
        nr_sed[i, k] = max(nr[i, k] + nr_tendency[i, k] * deltat, 0.0) / precip_fraction[i, k]
        qs_sed[i, k] = max(qs[i, k] + qs_tendency[i, k] * deltat, 0.0) / precip_fraction[i, k]
        ns_sed[i, k] = max(ns[i, k] + ns_tendency[i, k] * deltat, 0.0) / precip_fraction[i, k]
        # 🟤 graupel stays a grid mean here (FINDINGS.md)
        qg_sed[i, k] = max(qg[i, k] + qg_tendency[i, k] * deltat, 0.0)
        ng_sed[i, k] = max(ng[i, k] + ng_tendency[i, k] * deltat, 0.0)
        if settings.constant_droplet_number
            nc_sed[i, k] = settings.droplet_number / ρ[i, k]
        end
        if settings.constant_ice_number
            ni_sed[i, k] = settings.ice_number / ρ[i, k]
        end
        if settings.constant_graupel_number
            ng_sed[i, k] = settings.graupel_number / ρ[i, k] * precip_fraction[i, k]
        end
        qc_sed[i, k] = min(qc_sed[i, k], 5.0e-3)
        qi_sed[i, k] = min(qi_sed[i, k], 5.0e-3)
        qr_sed[i, k] = min(qr_sed[i, k], 10.0e-3)
        qs_sed[i, k] = min(qs_sed[i, k], 10.0e-3)
        qg_sed[i, k] = min(qg_sed[i, k], 10.0e-3)
    end

    # 🧊 cloud ice
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            number_before[i, k] = ni_sed[i, k]
            ice = size_distribution(hydrometeors.cloud_ice, qi_sed[i, k], ni_sed[i, k])
            λi[i, k], ni_sed[i, k], ni_intercept[i, k] = ice.slope, ice.number, ice.intercept
        end
        for k in 1:levels, i in 1:columns
            if qi_sed[i, k] >= qsmall
                if ni_sed[i, k] != number_before[i, k]
                    ni_tendency[i, k] = (ni_sed[i, k] * ice_fraction[i, k] - ni[i, k]) * rdeltat
                end
                radius_cloud_ice[i, k] = 1.5 / λi[i, k] * 1.0e6
                radius_cloud_ice[i, k] = radius_cloud_ice[i, k] * effi_factor
                x = λi[i, k]
                area_density_ice[i, k] = 2.0 * π * (1.0 / (x * x * x)) * ni_intercept[i, k] * ρ[i, k] *
                               1.0e-2
            else
                radius_cloud_ice[i, k] = 25.0
                radius_cloud_ice[i, k] = radius_cloud_ice[i, k] * effi_factor
                area_density_ice[i, k] = 0.0
            end
            diameter_cloud_ice[i, k] = radius_cloud_ice[i, k] * rhoi / rhows * 2.0
        end
    else
        for k in 1:levels, i in 1:columns
            radius_cloud_ice[i, k] = ice_radius_in[i, k] * 1.0e6
            radius_cloud_ice[i, k] = radius_cloud_ice[i, k] * effi_factor
            diameter_cloud_ice[i, k] = radius_cloud_ice[i, k] * 2.0
            area_density_ice[i, k] = 4.0 * π * (radius_cloud_ice[i, k] * radius_cloud_ice[i, k]) * ni[i, k] * ρ[i, k] * 1e-2
        end
    end

    # 💧 droplets, then again with a fixed number of 10⁸ kg⁻¹
    for k in 1:levels, i in 1:columns
        number_before[i, k] = nc_sed[i, k]
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, qc_sed[i, k],
            nc_sed[i, k], ρ[i, k])
        μc[i, k], λc[i, k], nc_sed[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        if qc_sed[i, k] >= qsmall
            if settings.constant_droplet_number
                nc_tendency[i, k] = (settings.droplet_number / ρ[i, k] * liquid_fraction[i, k] - nc[i, k]) *
                               rdeltat
            end
            if number_before[i, k] != nc_sed[i, k]
                nc_tendency[i, k] = (nc_sed[i, k] * liquid_fraction[i, k] - nc[i, k]) * rdeltat
            end
            radius_cloud_liquid[i, k] = (μc[i, k] + 3.0) / λc[i, k] / 2.0 * 1.0e6
            λc_radiation[i, k] = λc[i, k]
            μc_radiation[i, k] = μc[i, k]
            nc_sed[i, k] = 1.0e8
        end
    end
    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, qc_sed[i, k],
            nc_sed[i, k], ρ[i, k])
        μc[i, k], λc[i, k], nc_sed[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        if qc_sed[i, k] >= qsmall
            radius_cloud_liquid_fixed[i, k] = (μc[i, k] + 3.0) / λc[i, k] / 2.0 * 1.0e6
        else
            radius_cloud_liquid[i, k] = 10.0
            λc_radiation[i, k] = 0.0
            μc_radiation[i, k] = 0.0
            radius_cloud_liquid_fixed[i, k] = 10.0
        end
        number_before[i, k] = nr_sed[i, k]
    end

    # 🌧️ rain
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qr_sed[i, k], nr_sed[i, k])
        λr[i, k], nr_sed[i, k], N₀r[i, k] = rain.slope, rain.number, rain.intercept
    end
    for k in 1:levels, i in 1:columns
        if qr_sed[i, k] >= qsmall && number_before[i, k] != nr_sed[i, k]
            nr_tendency[i, k] = (nr_sed[i, k] * precip_fraction[i, k] - nr[i, k]) * rdeltat
        end
        number_before[i, k] = ns_sed[i, k]
    end

    # ❄️ snow
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, qs_sed[i, k], ns_sed[i, k])
        λs[i, k], ns_sed[i, k], ns_intercept[i, k] = snow.slope, snow.number, snow.intercept
    end
    for k in 1:levels, i in 1:columns
        if qs_sed[i, k] >= qsmall
            if number_before[i, k] != ns_sed[i, k]
                ns_tendency[i, k] = (ns_sed[i, k] * precip_fraction[i, k] - ns[i, k]) * rdeltat
            end
            x = λs[i, k]
            area_density_snow[i, k] = 2.0 * π * (1.0 / (x * x * x)) * ns_intercept[i, k] * ρ[i, k] *
                            1.0e-2
        end
        number_before[i, k] = ng_sed[i, k]
    end

    # 🟤 graupel or hail (hail first, then graupel, as the Fortran calls them)
    for (on, props) in ((settings.do_hail, hydrometeors.hail),
        (settings.do_graupel, hydrometeors.graupel))
        on || continue
        for k in 1:levels, i in 1:columns
            dense = size_distribution(props, qg_sed[i, k], ng_sed[i, k])
            λg[i, k], ng_sed[i, k] = dense.slope, dense.number
        end
    end

    for k in 1:levels, i in 1:columns
        if qg_sed[i, k] >= qsmall && number_before[i, k] != ng_sed[i, k]
            ng_tendency[i, k] = (ng_sed[i, k] * precip_fraction[i, k] - ng[i, k]) * rdeltat
        end
        # 🚫 no number without mass
        if qc[i, k] + qc_tendency[i, k] * deltat < qsmall
            nc_tendency[i, k] = -nc[i, k] * rdeltat
        end
        if settings.prognostic_cloud_ice && qi[i, k] + qi_tendency[i, k] * deltat < qsmall
            ni_tendency[i, k] = -ni[i, k] * rdeltat
        end
        if qr[i, k] + qr_tendency[i, k] * deltat < qsmall
            nr_tendency[i, k] = -nr[i, k] * rdeltat
        end
        if qs[i, k] + qs_tendency[i, k] * deltat < qsmall
            ns_tendency[i, k] = -ns[i, k] * rdeltat
        end
        if qg[i, k] + qg_tendency[i, k] * deltat < qsmall
            ng_tendency[i, k] = -ng[i, k] * rdeltat
        end
        # 📤 qc and qi are only used for output from here on
        qc[i, k] = qc[i, k] + qc_tendency[i, k] * deltat
        qi[i, k] = qi[i, k] + qi_tendency[i, k] * deltat
        if nr_per_volume[i, k] == 0.0
            nr_per_volume[i, k] = 1.0e-34
        end
    end
    return nothing
end
