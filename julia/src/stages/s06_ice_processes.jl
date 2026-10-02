# ❄️🧊 Stage 6, process rates of ice, snow and mixed phase (micro_pumas_v1.F90 lines
#    2166–2530).

"Scalar air NamedTuple of the vapor diffusion twins at point (i, k)."
function air_at(work, i, k)
    (; inputs, state) = work
    (; air) = state
    return (T = inputs.T[i, k], ρ = air.density[i, k], diffusivity = air.diffusivity[i, k],
        viscosity = air.viscosity[i, k], schmidt = air.schmidt[i, k],
        q_vapor = inputs.q_vapor[i, k], q_sat_liquid = air.q_sat_liquid[i, k],
        q_sat_ice = air.q_sat_ice[i, k])
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
    (; inputs, state, incloud, rates, process_rates, scratch) = work
    qs_from_host = inputs.snow_from_host.q  # tnd_qsnow
    ns_from_host = inputs.snow_from_host.N  # tnd_nsnow
    immersion_freezing_host = inputs.host_freezing.immersion  # frzimm
    contact_freezing_host = inputs.host_freezing.contact  # frzcnt
    deposition_freezing_host = inputs.host_freezing.deposition  # frzdep
    T = inputs.T  # t
    qv = inputs.q_vapor  # q
    ρ = state.air.density  # rho
    η = state.air.viscosity  # mu
    ρ_correction = state.air.density_correction  # rhof
    a_rain = state.air.fall_speed_prefactor.rain  # arn
    a_snow = state.air.fall_speed_prefactor.snow  # asn
    a_graupel = state.air.fall_speed_prefactor.graupel  # agn
    relative_variance = inputs.relative_variance  # relvar
    accretion_enhancement = inputs.accretion_enhancement  # accre_enhan
    cloud_fraction = incloud.cloud_fraction  # cldm
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    qi = state.cloud_ice.q
    ni = state.cloud_ice.N
    qs = state.snow.q
    ns = state.snow.N
    qg = state.graupel.q
    ng = state.graupel.N
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    λ_power = scratch.number_before  # dum_2D, here λ to the fall speed exponent
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    μc = incloud.cloud_liquid.μ  # pgam
    λc = incloud.cloud_liquid.λ  # lamc
    qi_cloud = incloud.cloud_ice.q  # qiic
    ni_cloud = incloud.cloud_ice.N  # niic
    λi = incloud.cloud_ice.λ  # lami
    N₀i = incloud.cloud_ice.N₀  # n0i
    qr_precip = incloud.rain.q  # qric
    nr_precip = incloud.rain.N  # nric
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    vr_number = incloud.rain.fall_speed_number  # unr
    qs_precip = incloud.snow.q  # qsic
    ns_precip = incloud.snow.N  # nsic
    λs = incloud.snow.λ  # lams
    N₀s = incloud.snow.N₀  # n0s
    vs_number = incloud.snow.fall_speed_number  # uns
    qg_precip = incloud.graupel.q  # qgic
    ng_precip = incloud.graupel.N  # ngic
    λg = incloud.graupel.λ  # lamg
    N₀g = incloud.graupel.N₀  # n0g
    vg_number = incloud.graupel.fall_speed_number  # ung
    vr_mass = process_rates.rain.fall_speed.mass  # umr
    vs_mass = process_rates.snow.fall_speed.mass  # ums
    vg_mass = process_rates.graupel.fall_speed.mass  # umg
    qi_autoconversion = rates.cloud_ice.autoconversion_to_snow.mass  # prci
    ni_autoconversion = rates.cloud_ice.autoconversion_to_snow.number  # nprci
    qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc
    nr_autoconversion = rates.rain.autoconversion_from_cloud_liquid.number  # nprc
    qc_accretion = rates.cloud_liquid.accretion_by_rain.mass  # pra
    nc_accretion = rates.cloud_liquid.accretion_by_rain.number  # npra
    nr_self_collection = rates.rain.self_collection.number  # nragg
    ns_aggregation = rates.snow.aggregation.number  # nsagg
    qc_riming_snow = rates.cloud_liquid.riming_by_snow.mass  # psacws
    nc_riming_snow = rates.cloud_liquid.riming_by_snow.number  # npsacws
    qi_splintering_snow = rates.cloud_ice.splintering_by_snow.mass  # msacwi
    ni_splintering_snow = rates.cloud_ice.splintering_by_snow.number  # nsacwi
    qr_accretion_snow = rates.rain.accretion_by_snow.mass  # pracs
    nr_accretion_snow = rates.rain.accretion_by_snow.number  # npracs
    qr_freezing_snow = rates.rain.freezing_to_snow.mass  # mnuccr
    nr_freezing_snow = rates.rain.freezing_to_snow.number  # nnuccr
    qc_immersion_freezing = rates.cloud_liquid.immersion_freezing.mass  # mnuccc
    nc_immersion_freezing = rates.cloud_liquid.immersion_freezing.number  # nnuccc
    qc_contact_freezing = rates.cloud_liquid.contact_freezing.mass  # mnucct
    nc_contact_freezing = rates.cloud_liquid.contact_freezing.number  # nnucct
    qc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.mass  # mnudep
    nc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.number  # nnudep
    qi_accretion_snow = rates.cloud_ice.accretion_by_snow.mass  # prai
    ni_accretion_snow = rates.cloud_ice.accretion_by_snow.number  # nprai
    qc_bergeron_snow = rates.cloud_liquid.bergeron_to_snow.mass  # bergs
    qs_deposition = rates.snow.deposition.mass  # vap_deps
    qi_deposition = rates.cloud_ice.deposition.mass  # vap_dep
    qi_sublimation = rates.cloud_ice.sublimation.mass  # ice_sublim
    qc_bergeron_ice = rates.cloud_liquid.bergeron_to_cloud_ice.mass  # berg
    ni_sublimation = rates.cloud_ice.sublimation.number  # nsubi
    nc_evaporation = rates.cloud_liquid.evaporation.number  # nsubc
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        ice = size_distribution(hydrometeors.cloud_ice, qi_cloud[i, k], ni_cloud[i, k])
        λi[i, k], N₀i[i, k], ni_cloud[i, k] = ice.slope, ice.intercept, ice.number
    end

    # 🧊➡️❄️ autoconversion of cloud ice to snow
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            rate = ice_to_snow_autoconversion(
                parameters, T[i, k], qi_cloud[i, k], λi[i, k],
                N₀i[i, k])
            qi_autoconversion[i, k], ni_autoconversion[i, k] = rate.mass, rate.number
        else
            qi_autoconversion[i, k] = qs_from_host[i, k] / cloud_fraction[i, k]
            ni_autoconversion[i, k] = ns_from_host[i, k] / cloud_fraction[i, k]
        end
    end

    # ❄️🟤 in-precipitation snow and graupel
    for k in 1:levels, i in 1:columns
        qs_precip[i, k] = qs[i, k] / precip_fraction[i, k]
        ns_precip[i, k] = ns[i, k] / precip_fraction[i, k]
        qs_precip[i, k] = min(qs_precip[i, k], 0.01)
        if qs_precip[i, k] < smallest_mixing_ratio
            qs_precip[i, k] = 0.0
            ns_precip[i, k] = 0.0
        end
        ns_precip[i, k] = max(ns_precip[i, k], 0.0)
        qg_precip[i, k] = qg[i, k] / precip_fraction[i, k]
        ng_precip[i, k] = ng[i, k] / precip_fraction[i, k]
        qg_precip[i, k] = min(qg_precip[i, k], 0.01)
        if qg_precip[i, k] < smallest_mixing_ratio
            qg_precip[i, k] = 0.0
            ng_precip[i, k] = 0.0
        end
        ng_precip[i, k] = max(ng_precip[i, k], 0.0)
    end

    # 🌧️ provisional rain fall speeds (m/s)
    br = fall_speed_law.rain.b
    for k in 1:levels, i in 1:columns
        if λr[i, k] >= smallest_mixing_ratio
            λ_power[i, k] = λr[i, k]^br
            vr_number[i, k] = min(
                a_rain[i, k] * gamma.rain_plus1 / λ_power[i, k], 9.1 * ρ_correction[i, k])
            vr_mass[i, k] = min(a_rain[i, k] * gamma.rain_plus4 / (6.0 * λ_power[i, k]),
                9.1 * ρ_correction[i, k])
        else
            vr_mass[i, k] = 0.0
            vr_number[i, k] = 0.0
        end
    end

    # ❄️ snow size distribution and provisional fall speeds
    for k in 1:levels, i in 1:columns
        snow = size_distribution(hydrometeors.snow, qs_precip[i, k], ns_precip[i, k])
        λs[i, k], N₀s[i, k], ns_precip[i, k] = snow.slope, snow.intercept, snow.number
    end
    bs = fall_speed_law.snow.b
    for k in 1:levels, i in 1:columns
        if settings.ifs_snow_fall_speed
            vs_mass[i, k] = λs[i, k] > 0.0 ? 1.0 : 0.0
            vs_number[i, k] = λs[i, k] > 0.0 ? 1.0 : 0.0
        elseif λs[i, k] > 0.0
            λ_power[i, k] = λs[i, k]^bs
            vs_mass[i, k] = min(a_snow[i, k] * gamma.snow_plus4 / (6.0 * λ_power[i, k]),
                1.2 * ρ_correction[i, k])
            vs_mass[i, k] = vs_mass[i, k] * settings.snow_fall_speed_factor
            vs_number[i, k] = min(
                a_snow[i, k] * gamma.snow_plus1 / λ_power[i, k], 1.2 * ρ_correction[i, k])
        else
            vs_mass[i, k] = 0.0
            vs_number[i, k] = 0.0
        end
    end

    # 🟤 graupel (or hail) size distribution and provisional fall speeds
    if settings.do_hail || settings.do_graupel
        dense_props = settings.do_hail ? hydrometeors.hail : hydrometeors.graupel
        for k in 1:levels, i in 1:columns
            graupel = size_distribution(dense_props, qg_precip[i, k], ng_precip[i, k])
            λg[i, k], N₀g[i, k], ng_precip[i, k] = graupel.slope, graupel.intercept,
            graupel.number
        end
    end
    for k in 1:levels, i in 1:columns
        if λg[i, k] > 0.0
            λ_power[i, k] = λg[i, k]^dense.b
            vg_mass[i, k] = min(
                a_graupel[i, k] * gamma.graupel_plus4 / (6.0 * λ_power[i, k]),
                20.0 * ρ_correction[i, k])
            vg_number[i, k] = min(a_graupel[i, k] * gamma.graupel_plus1 / λ_power[i, k],
                20.0 * ρ_correction[i, k])
        else
            vg_mass[i, k] = 0.0
            vg_number[i, k] = 0.0
        end
    end

    # 💧➡️🧊 droplet freezing from the host's heterogeneous freezing rates
    if settings.prognostic_cloud_ice
        settings.host_heterogeneous_freezing ||
            error("Bigg droplet freezing (hetfrz off) is not translated yet")
        for k in 1:levels, i in 1:columns
            mi0l = qc_cloud[i, k] / max(nc_cloud[i, k], 1.0e6 / ρ[i, k])
            mi0l = max(minimum_frozen_droplet_mass, mi0l)
            if qc_cloud[i, k] >= smallest_mixing_ratio
                nc_immersion_freezing[i, k] = immersion_freezing_host[i, k] * 1.0e6 /
                                              ρ[i, k]
                qc_immersion_freezing[i, k] = nc_immersion_freezing[i, k] * mi0l
                nc_contact_freezing[i, k] = contact_freezing_host[i, k] * 1.0e6 / ρ[i, k]
                qc_contact_freezing[i, k] = nc_contact_freezing[i, k] * mi0l
                nc_deposition_nucleation[i, k] = deposition_freezing_host[i, k] * 1.0e6 /
                                                 ρ[i, k]
                qc_deposition_nucleation[i, k] = nc_deposition_nucleation[i, k] *
                                                 new_ice_crystal_mass
            else
                nc_immersion_freezing[i, k] = 0.0
                qc_immersion_freezing[i, k] = 0.0
                nc_contact_freezing[i, k] = 0.0
                qc_contact_freezing[i, k] = 0.0
                nc_deposition_nucleation[i, k] = 0.0
                qc_deposition_nucleation[i, k] = 0.0
            end
        end
    else
        for array in (qc_immersion_freezing, nc_immersion_freezing, qc_contact_freezing,
            nc_contact_freezing, qc_deposition_nucleation, nc_deposition_nucleation)
            fill!(array, 0.0)
        end
    end

    iaccr = settings.ice_accretion_factor
    for k in 1:levels, i in 1:columns
        ns_aggregation[i, k] = snow_self_aggregation(
            parameters, T[i, k], ρ[i, k], a_snow[i, k],
            ρ_snow, qs_precip[i, k], ns_precip[i, k])
    end
    for k in 1:levels, i in 1:columns
        riming = droplet_collection_by_snow(parameters, T[i, k], ρ[i, k], a_snow[i, k],
            vs_number[i, k], η[i, k], qc_cloud[i, k], nc_cloud[i, k], qs_precip[i, k], μc[
                i, k],
            λc[i, k], λs[i, k], N₀s[i, k])
        qc_riming_snow[i, k], nc_riming_snow[i, k] = riming.mass, riming.number
    end
    for k in 1:levels, i in 1:columns
        qc_riming_snow[i, k] = qc_riming_snow[i, k] * iaccr
        nc_riming_snow[i, k] = nc_riming_snow[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            splinters = rime_splintering_by_snow(T[i, k], qc_riming_snow[i, k])
            qi_splintering_snow[i, k], ni_splintering_snow[i, k] = splinters.mass,
            splinters.number
            qc_riming_snow[i, k] = splinters.riming
        else
            ni_splintering_snow[i, k] = 0.0
            qi_splintering_snow[i, k] = 0.0
        end
    end
    for k in 1:levels, i in 1:columns
        collection = rain_collection_by_snow(parameters, T[i, k], ρ[i, k], vr_mass[i, k],
            vs_mass[i, k], vr_number[i, k], vs_number[i, k], qr_precip[i, k],
            qs_precip[i, k], λr[i, k], N₀r[i, k],
            λs[i, k], N₀s[i, k])
        qr_accretion_snow[i, k], nr_accretion_snow[i, k] = collection.mass,
        collection.number
    end
    for k in 1:levels, i in 1:columns
        freezing = rain_immersion_freezing(
            parameters, T[i, k], qr_precip[i, k], nr_precip[i, k],
            λr[i, k])
        qr_freezing_snow[i, k], nr_freezing_snow[i, k] = freezing.mass, freezing.number
    end

    # 💧🌧️ accretion of cloud liquid by rain
    if settings.warm_rain == :sb2001
        qc_sb2001 = process_rates.cloud_liquid.warm_rain_sb2001.mass  # qctend_SB2001
        nc_sb2001 = process_rates.cloud_liquid.warm_rain_sb2001.number  # nctend_SB2001
        qr_sb2001 = process_rates.rain.warm_rain_sb2001.mass  # qrtend_SB2001
        nr_sb2001 = process_rates.rain.warm_rain_sb2001.number  # nrtend_SB2001
        for k in 1:levels, i in 1:columns
            rate = sb2001_accretion(
                qc_cloud[i, k], nc_cloud[i, k], qr_precip[i, k], ρ[i, k])
            qc_accretion[i, k], nc_accretion[i, k] = rate.mass, rate.droplet_number
        end
        for k in 1:levels, i in 1:columns
            nc_sb2001[i, k] = nc_sb2001[i, k] + nc_accretion[i, k]
            qc_sb2001[i, k] = qc_sb2001[i, k] + qc_accretion[i, k]
            nr_sb2001[i, k] = nr_sb2001[i, k] + nc_accretion[i, k]
            qr_sb2001[i, k] = qr_sb2001[i, k] - qc_accretion[i, k]
        end
    end
    if settings.warm_rain == :kk2000
        qc_kk2000 = process_rates.cloud_liquid.warm_rain_kk2000.mass  # qctend_KK2000
        nc_kk2000 = process_rates.cloud_liquid.warm_rain_kk2000.number  # nctend_KK2000
        qr_kk2000 = process_rates.rain.warm_rain_kk2000.mass  # qrtend_KK2000
        nr_kk2000 = process_rates.rain.warm_rain_kk2000.number  # nrtend_KK2000
        for k in 1:levels, i in 1:columns
            qr_for_accretion = qr_precip[i, k]
            qc_for_accretion = qc_cloud[i, k]
            nc_for_accretion = nc_cloud[i, k]
            if settings.accretion_sees_autoconversion
                qr_for_accretion = qr_for_accretion + qc_autoconversion[i, k] * Δt
                qc_for_accretion = qc_for_accretion - qc_autoconversion[i, k] * Δt
                nc_for_accretion = nc_for_accretion - nr_autoconversion[i, k] * Δt
            end
            rate = accretion(
                parameters, qr_for_accretion, qc_for_accretion, nc_for_accretion,
                relative_variance[i, k], accretion_enhancement[i, k])
            qc_accretion[i, k], nc_accretion[i, k] = rate.mass, rate.droplet_number
        end
        for k in 1:levels, i in 1:columns
            nc_kk2000[i, k] = nc_kk2000[i, k] + nc_accretion[i, k]
            qc_kk2000[i, k] = qc_kk2000[i, k] + qc_accretion[i, k]
            nr_kk2000[i, k] = nr_kk2000[i, k] + nc_accretion[i, k]
            qr_kk2000[i, k] = qr_kk2000[i, k] - qc_accretion[i, k]
        end
    end
    for k in 1:levels, i in 1:columns
        qc_accretion[i, k] = qc_accretion[i, k] * settings.accretion_enhancement
        nc_accretion[i, k] = nc_accretion[i, k] * settings.accretion_enhancement
    end
    for k in 1:levels, i in 1:columns
        nr_self_collection[i, k] = rain_self_collection(
            ρ[i, k], qr_precip[i, k], nr_precip[i, k])
    end

    # 🧊❄️ accretion of cloud ice by snow
    for k in 1:levels, i in 1:columns
        if settings.prognostic_cloud_ice
            rate = ice_collection_by_snow(parameters, T[i, k], ρ[i, k], a_snow[i, k],
                qi_cloud[i, k], ni_cloud[i, k], qs_precip[i, k], λs[i, k], N₀s[i, k])
            qi_accretion_snow[i, k], ni_accretion_snow[i, k] = rate.mass, rate.number
        else
            qi_accretion_snow[i, k] = 0.0
            ni_accretion_snow[i, k] = 0.0
        end
    end

    # 🌫️ Bergeron process and vapor deposition onto snow and ice
    berg_eff = settings.bergeron_efficiency
    for k in 1:levels, i in 1:columns
        qc_bergeron_snow[i, k] = bergeron_to_snow(
            parameters, air_at(work, i, k), a_snow[i, k], qc_cloud[i, k],
            qs_precip[i, k], λs[i, k], N₀s[i, k])
    end
    for k in 1:levels, i in 1:columns
        qc_bergeron_snow[i, k] = qc_bergeron_snow[i, k] * berg_eff
    end
    for k in 1:levels, i in 1:columns
        qs_deposition[i, k] = snow_deposition(
            parameters, air_at(work, i, k), qs[i, k], ns[i, k],
            precip_fraction[i, k], a_snow[i, k])
    end
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            rate = ice_deposition(parameters, air_at(work, i, k), qi[i, k], ni[i, k],
                ice_fraction[i, k])
            qc_bergeron_ice[i, k], qi_deposition[i, k] = rate.bergeron, rate.deposition
            qi_sublimation[i, k] = rate.sublimation
        end
        for k in 1:levels, i in 1:columns
            qc_bergeron_ice[i, k] = qc_bergeron_ice[i, k] * berg_eff
            if qi_sublimation[i, k] < 0.0 && qi[i, k] > smallest_mixing_ratio &&
               ice_fraction[i, k] > minimum_cloud_fraction
                ni_sublimation[i, k] = number_sublimation_factor * qi_sublimation[i, k] /
                                       qi[i, k] * ni[i, k] /
                                       ice_fraction[i, k]
            else
                ni_sublimation[i, k] = 0.0
            end
            nc_evaporation[i, k] = 0.0
        end
    end
    return nothing
end
