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
    (; g, Rd) = host_constants(parameters)
    dense = dense_precipitation(parameters)
    (; inputs, state, incloud, process_rates, tendencies, scratch) = work
    T = inputs.T  # t
    p_interface = inputs.interface_pressure  # pint
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
    ρ_correction = state.air.density_correction  # rhof
    a_cloud_liquid = state.air.fall_speed_prefactor.cloud_liquid  # acn
    a_cloud_ice = state.air.fall_speed_prefactor.cloud_ice  # ain
    a_small_ice = state.air.fall_speed_prefactor.small_cloud_ice  # ajn
    a_rain = state.air.fall_speed_prefactor.rain  # arn
    a_snow = state.air.fall_speed_prefactor.snow  # asn
    a_graupel = state.air.fall_speed_prefactor.graupel  # agn
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
    ng_tendency = tendencies.graupel.N  # ngtend
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
    fall_qc = scratch.pressure_fall_speed.cloud_liquid.q  # fc
    fall_nc = scratch.pressure_fall_speed.cloud_liquid.N  # fnc
    fall_qi = scratch.pressure_fall_speed.cloud_ice.q  # fi
    fall_ni = scratch.pressure_fall_speed.cloud_ice.N  # fni
    fall_qr = scratch.pressure_fall_speed.rain.q  # fr
    fall_nr = scratch.pressure_fall_speed.rain.N  # fnr
    fall_qs = scratch.pressure_fall_speed.snow.q  # fs
    fall_ns = scratch.pressure_fall_speed.snow.N  # fns
    fall_qg = scratch.pressure_fall_speed.graupel.q  # fg
    fall_ng = scratch.pressure_fall_speed.graupel.N  # fng
    λi = incloud.cloud_ice.λ  # lami
    μc = incloud.cloud_liquid.μ  # pgam
    λc = incloud.cloud_liquid.λ  # lamc
    λr = incloud.rain.λ  # lamr
    λs = incloud.snow.λ  # lams
    λg = incloud.graupel.λ  # lamg
    vr_number = incloud.rain.fall_speed_number  # unr
    vs_number = incloud.snow.fall_speed_number  # uns
    vg_number = incloud.graupel.fall_speed_number  # ung
    z_interface = scratch.interface_height  # zint
    vc_mass = process_rates.cloud_liquid.fall_speed.mass  # vtrmc
    vi_mass = process_rates.cloud_ice.fall_speed.mass  # vtrmi
    vr_mass = process_rates.rain.fall_speed.mass  # umr
    vs_mass = process_rates.snow.fall_speed.mass  # ums
    vg_mass = process_rates.graupel.fall_speed.mass  # umg
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        qc_sed[i, k] = (qc[i, k] + qc_tendency[i, k] * Δt) / liquid_fraction[i, k]
        qi_sed[i, k] = (qi[i, k] + qi_tendency[i, k] * Δt) / ice_fraction[i, k]
        nc_sed[i, k] = max((nc[i, k] + nc_tendency[i, k] * Δt) / liquid_fraction[i, k], 0.0)
        ni_sed[i, k] = max((ni[i, k] + ni_tendency[i, k] * Δt) / ice_fraction[i, k], 0.0)
        qr_sed[i, k] = (qr[i, k] + qr_tendency[i, k] * Δt) / precip_fraction[i, k]
        nr_sed[i, k] = max((nr[i, k] + nr_tendency[i, k] * Δt) / precip_fraction[i, k], 0.0)
        qs_sed[i, k] = (qs[i, k] + qs_tendency[i, k] * Δt) / precip_fraction[i, k]
        ns_sed[i, k] = max((ns[i, k] + ns_tendency[i, k] * Δt) / precip_fraction[i, k], 0.0)
        qg_sed[i, k] = (qg[i, k] + qg_tendency[i, k] * Δt) / precip_fraction[i, k]
        ng_sed[i, k] = max((ng[i, k] + ng_tendency[i, k] * Δt) / precip_fraction[i, k], 0.0)
        # 🔢 constant numbers; the constant rain and snow numbers (nrcons, nscons) are
        #    not settings of oPUMAS.jl yet, so their branches are left out.
        if settings.constant_graupel_number
            ng_sed[i, k] = settings.graupel_number / ρ[i, k]
        end
        if settings.constant_droplet_number
            nc_sed[i, k] = settings.droplet_number / ρ[i, k]
        end
        if settings.constant_ice_number
            ni_sed[i, k] = settings.ice_number / ρ[i, k]
        end
    end

    for k in 1:levels, i in 1:columns
        ice = size_distribution(hydrometeors.cloud_ice, qi_sed[i, k], ni_sed[i, k])
        λi[i, k], ni_sed[i, k] = ice.slope, ice.number
    end
    for k in 1:levels, i in 1:columns
        liquid = droplet_size_distribution(hydrometeors.cloud_liquid, qc_sed[i, k],
            nc_sed[i, k], ρ[i, k])
        μc[i, k], λc[i, k], nc_sed[i, k] = liquid.shape, liquid.slope, liquid.number
    end
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qr_sed[i, k], nr_sed[i, k])
        λr[i, k], nr_sed[i, k] = rain.slope, rain.number
    end
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, qs_sed[i, k], ns_sed[i, k])
        λs[i, k], ns_sed[i, k] = snow.slope, snow.number
    end
    if settings.do_hail || settings.do_graupel
        dense_props = settings.do_hail ? hydrometeors.hail : hydrometeors.graupel
        for k in 1:levels, i in 1:columns
            graupel = size_distribution(dense_props, qg_sed[i, k], ng_sed[i, k])
            λg[i, k], ng_sed[i, k] = graupel.slope, graupel.number
        end
    end

    if settings.implicit_sedimentation
        for i in 1:columns
            z_interface[i, levels+1] = 0.0
            for k in levels:-1:1
                H = Rd * T[i, k] / g * log(p_interface[i, k+1] / p_interface[i, k])
                z_interface[i, k] = z_interface[i, k+1] + H
            end
        end
    end

    # 💧 droplets
    bc = fall_speed_law.cloud_liquid.b
    for k in 1:levels, i in 1:columns
        if qc_sed[i, k] >= smallest_mixing_ratio
            mass_gamma_top = 4.0 + bc + μc[i, k]
            mass_gamma_bottom = μc[i, k] + 4.0
            vc_mass[i, k] = a_cloud_liquid[i, k] * gamma_function(mass_gamma_top) /
                          (λc[i, k]^bc * gamma_function(mass_gamma_bottom))
            if settings.ifs_snow_fall_speed
                fall_qc[i, k] = 0.0
                fall_nc[i, k] = 0.0
            else
                number_gamma_top = 1.0 + bc + μc[i, k]
                number_gamma_bottom = μc[i, k] + 1.0
                fall_qc[i, k] = g * ρ[i, k] * vc_mass[i, k]
                fall_nc[i, k] = g * ρ[i, k] * a_cloud_liquid[i, k] * gamma_function(number_gamma_top) /
                            (λc[i, k]^bc * gamma_function(number_gamma_bottom))
            end
        else
            fall_qc[i, k] = 0.0
            fall_nc[i, k] = 0.0
        end
        qc_sed[i, k] = (qc[i, k] + qc_tendency[i, k] * Δt)
        nc_sed[i, k] = max((nc[i, k] + nc_tendency[i, k] * Δt), 0.0)
        if qc_sed[i, k] < smallest_mixing_ratio
            nc_sed[i, k] = 0.0
        end
    end

    # 🧊 cloud ice, blended with small-ice fall speeds below 20 μm
    bi, bj = fall_speed_law.cloud_ice.b, fall_speed_law.small_cloud_ice.b
    vtrmi_factor = settings.ice_fall_speed_factor
    for k in 1:levels, i in 1:columns
        if qi_sed[i, k] >= smallest_mixing_ratio
            vi_mass[i, k] = min(a_cloud_ice[i, k] * gamma.ice_plus4 / (6.0 * λi[i, k]^bi),
                1.2 * ρ_correction[i, k])
            vi_mass[i, k] = vi_mass[i, k] * vtrmi_factor
            fall_qi[i, k] = g * ρ[i, k] * vi_mass[i, k]
            fall_ni[i, k] = g * ρ[i, k] *
                        min(a_cloud_ice[i, k] * gamma.ice_plus1 / λi[i, k]^bi, 1.2 * ρ_correction[i, k])
            irad = 1.5 / λi[i, k] * 1e6
            ifrac = min(1.0, max(0.0, (irad - 18.0) / 2.0))
            if ifrac < 1.0
                # 🔁 the factor is applied again to the blend (FINDINGS.md)
                vi_mass[i, k] = ifrac * vi_mass[i, k] +
                              (1.0 - ifrac) *
                              min(a_small_ice[i, k] * gamma.small_ice_plus4 / (6.0 * λi[i, k]^bj),
                    1.2 * ρ_correction[i, k])
                vi_mass[i, k] = vi_mass[i, k] * vtrmi_factor
                fall_qi[i, k] = g * ρ[i, k] * vi_mass[i, k]
                fall_ni[i, k] = ifrac * fall_ni[i, k] +
                            (1.0 - ifrac) * g * ρ[i, k] *
                            min(a_small_ice[i, k] * gamma.small_ice_plus1 / λi[i, k]^bj,
                    1.2 * ρ_correction[i, k])
            end
            if settings.ifs_snow_fall_speed
                fall_qi[i, k] = g * ρ[i, k] * 0.1
                fall_ni[i, k] = g * ρ[i, k] * 0.1
            end
        else
            fall_qi[i, k] = 0.0
            fall_ni[i, k] = 0.0
        end
        qi_sed[i, k] = (qi[i, k] + qi_tendency[i, k] * Δt)
        ni_sed[i, k] = max((ni[i, k] + ni_tendency[i, k] * Δt), 0.0)
        if qi_sed[i, k] < smallest_mixing_ratio
            ni_sed[i, k] = 0.0
        end
    end

    # 🌧️ rain; the fall correction starts at k = 3 as in the Fortran (FINDINGS.md)
    br = fall_speed_law.rain.b
    for i in 1:columns, k in 1:levels
        if λr[i, k] >= smallest_mixing_ratio
            λ_power = λr[i, k]^br
            vr_number[i, k] = min(a_rain[i, k] * gamma.rain_plus1 / λ_power, 9.1 * ρ_correction[i, k])
            fall_nr[i, k] = g * ρ[i, k] * vr_number[i, k]
            vr_mass[i, k] = min(a_rain[i, k] * gamma.rain_plus4 / (6.0 * λ_power), 9.1 * ρ_correction[i, k])
            fall_qr[i, k] = g * ρ[i, k] * vr_mass[i, k]
        else
            fall_qr[i, k] = 0.0
            fall_nr[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fall_qr[i, k] < 1.e-10
            fall_qr[i, k] = fall_qr[i, k-1]
            fall_nr[i, k] = fall_nr[i, k-1]
        end
        qr_sed[i, k] = (qr[i, k] + qr_tendency[i, k] * Δt)
        nr_sed[i, k] = max((nr[i, k] + nr_tendency[i, k] * Δt), 0.0)
        if qr_sed[i, k] < smallest_mixing_ratio
            nr_sed[i, k] = 0.0
        end
    end

    # ❄️ snow
    bs = fall_speed_law.snow.b
    for i in 1:columns, k in 1:levels
        if λs[i, k] >= smallest_mixing_ratio
            λ_power = λs[i, k]^bs
            vs_mass[i, k] = min(a_snow[i, k] * gamma.snow_plus4 / (6.0 * λ_power), 1.2 * ρ_correction[i, k])
            vs_mass[i, k] = vs_mass[i, k] * settings.snow_fall_speed_factor
            fall_qs[i, k] = g * ρ[i, k] * vs_mass[i, k]
            vs_number[i, k] = min(a_snow[i, k] * gamma.snow_plus1 / λ_power, 1.2 * ρ_correction[i, k])
            fall_ns[i, k] = g * ρ[i, k] * vs_number[i, k]
            if settings.ifs_snow_fall_speed
                vs_mass[i, k] = 1.0
                vs_number[i, k] = 1.0
            end
        else
            fall_qs[i, k] = 0.0
            fall_ns[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fall_qs[i, k] < 1.e-10
            fall_qs[i, k] = fall_qs[i, k-1]
            fall_ns[i, k] = fall_ns[i, k-1]
        end
        qs_sed[i, k] = (qs[i, k] + qs_tendency[i, k] * Δt)
        ns_sed[i, k] = max((ns[i, k] + ns_tendency[i, k] * Δt), 0.0)
        if qs_sed[i, k] < smallest_mixing_ratio
            ns_sed[i, k] = 0.0
        end
    end

    # 🟤 graupel or hail
    for i in 1:columns, k in 1:levels
        if λg[i, k] >= smallest_mixing_ratio
            λ_power = λg[i, k]^dense.b
            vg_mass[i, k] = min(a_graupel[i, k] * gamma.graupel_plus4 / (6.0 * λ_power),
                20.0 * ρ_correction[i, k])
            fall_qg[i, k] = g * ρ[i, k] * vg_mass[i, k]
            vg_number[i, k] = min(a_graupel[i, k] * gamma.graupel_plus1 / λ_power, 20.0 * ρ_correction[i, k])
            fall_ng[i, k] = g * ρ[i, k] * vg_number[i, k]
        else
            fall_qg[i, k] = 0.0
            fall_ng[i, k] = 0.0
        end
        if settings.precipitation_fall_correction && k > 2 && fall_qg[i, k] < 1.e-10
            fall_qg[i, k] = fall_qg[i, k-1]
            fall_ng[i, k] = fall_ng[i, k-1]
        end
        qg_sed[i, k] = (qg[i, k] + qg_tendency[i, k] * Δt)
        ng_sed[i, k] = max((ng[i, k] + ng_tendency[i, k] * Δt), 0.0)
        if qg_sed[i, k] < smallest_mixing_ratio
            ng_sed[i, k] = 0.0
        end
    end
    return nothing
end
