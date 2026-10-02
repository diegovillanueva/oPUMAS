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
    (; cp, Rv, Lv, Ls) = host_constants(parameters)
    Lv_squared, Ls_squared = Lv^2, Ls^2
    (; inputs, state, incloud, rates, diagnostics, scratch) = work
    T = inputs.T  # t
    qv = inputs.q_vapor  # q
    p = inputs.pressure
    qsat_ice = state.air.q_sat_ice  # qvi
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    qr = state.rain.q
    nr = state.rain.N
    qi = state.cloud_ice.q
    ni = state.cloud_ice.N
    qs = state.snow.q
    ns = state.snow.N
    qg = state.graupel.q
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    λr = incloud.rain.λ  # lamr
    qc_limiter_ratio = diagnostics.cloud_liquid.limiter_ratio  # qcrat
    qv_after_deposition = scratch.number_before  # dum_2D, here the vapor after deposition
    T_estimate = scratch.temperature_estimate  # ttmpA
    qv_estimate = scratch.vapor_estimate  # qtmpAI
    qsat_estimate = scratch.q_sat_estimate  # qvnAI
    qsat_adjusted = scratch.q_sat_adjusted  # qvnA
    rain_share = scratch.evaporation_share.rain  # dum1A
    snow_share = scratch.evaporation_share.snow  # dum2A
    graupel_share = scratch.evaporation_share.graupel  # dum3A
    qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc
    nr_autoconversion = rates.rain.autoconversion_from_cloud_liquid.number  # nprc
    nc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.number  # nprc1
    qc_accretion = rates.cloud_liquid.accretion_by_rain.mass  # pra
    nc_accretion = rates.cloud_liquid.accretion_by_rain.number  # npra
    nr_self_collection = rates.rain.self_collection.number  # nragg
    ns_aggregation = rates.snow.aggregation.number  # nsagg
    qc_riming_snow = rates.cloud_liquid.riming_by_snow.mass  # psacws
    nc_riming_snow = rates.cloud_liquid.riming_by_snow.number  # npsacws
    qr_accretion_snow = rates.rain.accretion_by_snow.mass  # pracs
    nr_accretion_snow = rates.rain.accretion_by_snow.number  # npracs
    qc_immersion_freezing = rates.cloud_liquid.immersion_freezing.mass  # mnuccc
    nc_immersion_freezing = rates.cloud_liquid.immersion_freezing.number  # nnuccc
    qc_contact_freezing = rates.cloud_liquid.contact_freezing.mass  # mnucct
    nc_contact_freezing = rates.cloud_liquid.contact_freezing.number  # nnucct
    qc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.mass  # mnudep
    nc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.number  # nnudep
    qi_nucleation = rates.cloud_ice.nucleation_from_vapor.mass  # mnuccd
    ni_nucleation = rates.cloud_ice.nucleation_from_vapor.number  # nnuccd
    qi_splintering_snow = rates.cloud_ice.splintering_by_snow.mass  # msacwi
    ni_splintering_snow = rates.cloud_ice.splintering_by_snow.number  # nsacwi
    qr_freezing_snow = rates.rain.freezing_to_snow.mass  # mnuccr
    nr_freezing_snow = rates.rain.freezing_to_snow.number  # nnuccr
    qr_freezing_ice = rates.rain.freezing_to_cloud_ice.mass  # mnuccri
    nr_freezing_ice = rates.rain.freezing_to_cloud_ice.number  # nnuccri
    qi_autoconversion = rates.cloud_ice.autoconversion_to_snow.mass  # prci
    ni_autoconversion = rates.cloud_ice.autoconversion_to_snow.number  # nprci
    qi_accretion_snow = rates.cloud_ice.accretion_by_snow.mass  # prai
    ni_accretion_snow = rates.cloud_ice.accretion_by_snow.number  # nprai
    qr_evaporation = rates.rain.evaporation.mass  # pre
    nr_evaporation = rates.rain.evaporation.number  # nsubr
    qs_sublimation = rates.snow.sublimation.mass  # prds
    ns_sublimation = rates.snow.sublimation.number  # nsubs
    qi_deposition = rates.cloud_ice.deposition.mass  # vap_dep
    qi_sublimation = rates.cloud_ice.sublimation.mass  # ice_sublim
    ni_sublimation = rates.cloud_ice.sublimation.number  # nsubi
    qs_deposition = rates.snow.deposition.mass  # vap_deps
    nc_evaporation = rates.cloud_liquid.evaporation.number  # nsubc
    qc_bergeron_ice = rates.cloud_liquid.bergeron_to_cloud_ice.mass  # berg
    qc_bergeron_snow = rates.cloud_liquid.bergeron_to_snow.mass  # bergs
    qr_accretion_graupel = rates.rain.accretion_by_graupel.mass  # pracg
    nr_accretion_graupel = rates.rain.accretion_by_graupel.number  # npracg
    qc_riming_graupel = rates.cloud_liquid.riming_by_graupel.mass  # psacwg
    nc_riming_graupel = rates.cloud_liquid.riming_by_graupel.number  # npsacwg
    qg_rimed_snow = rates.graupel.from_snow_riming.mass  # pgsacw
    ng_rimed_snow = rates.graupel.from_snow_riming.number  # nscng
    qg_rain_snow = rates.graupel.rain_snow_collection.mass  # pgracs
    ng_rain_snow = rates.graupel.rain_snow_collection.number  # ngracs
    qs_collection_rain = rates.snow.collection_by_rain.mass  # psacr
    qg_deposition = rates.graupel.deposition.mass  # prdg
    qi_splintering_graupel = rates.cloud_ice.splintering_by_graupel.mass  # qmultg
    ni_splintering_graupel = rates.cloud_ice.splintering_by_graupel.number  # nmultg
    qi_splintering_rain = rates.cloud_ice.splintering_rain_graupel.mass  # qmultrg
    ni_splintering_rain = rates.cloud_ice.splintering_rain_graupel.number  # nmultrg
    Δt = work.inputs.time_step
    inv_Δt = 1.0 / Δt
    graupel_on = settings.do_hail || settings.do_graupel
    columns, levels = grid_size(work)

    if settings.ifs_evaporation_scaling
        for k in 1:levels, i in 1:columns
            qr_evaporation[i, k] = 0.15 * qr_evaporation[i, k]
        end
    end

    # 💧 conservation of qc, then the deposition limit
    for k in 1:levels, i in 1:columns
        sink_over_step = ((qc_autoconversion[i, k] + qc_accretion[i, k] +
                           qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                           qi_splintering_snow[i, k] +
                           qc_riming_snow[i, k] + qc_bergeron_snow[i, k] +
                           qi_splintering_graupel[i, k] + qc_riming_graupel[i, k] +
                           qg_rimed_snow[i, k]) *
                          liquid_fraction[i, k] + qc_bergeron_ice[i, k]) * Δt
        if sink_over_step > qc[i, k]
            ratio = qc[i, k] * inv_Δt /
                    ((qc_autoconversion[i, k] + qc_accretion[i, k] +
                      qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                      qi_splintering_snow[i, k] +
                      qc_riming_snow[i, k] + qc_bergeron_snow[i, k] +
                      qi_splintering_graupel[i, k] + qc_riming_graupel[i, k] +
                      qg_rimed_snow[i, k]) * liquid_fraction[i, k] +
                     qc_bergeron_ice[i, k]) * one_minus_small
            qi_splintering_graupel[i, k] = qi_splintering_graupel[i, k] * ratio
            qc_riming_graupel[i, k] = qc_riming_graupel[i, k] * ratio
            qg_rimed_snow[i, k] = qg_rimed_snow[i, k] * ratio
            qc_autoconversion[i, k] = qc_autoconversion[i, k] * ratio
            qc_accretion[i, k] = qc_accretion[i, k] * ratio
            qc_immersion_freezing[i, k] = qc_immersion_freezing[i, k] * ratio
            qc_contact_freezing[i, k] = qc_contact_freezing[i, k] * ratio
            qi_splintering_snow[i, k] = qi_splintering_snow[i, k] * ratio
            qc_riming_snow[i, k] = qc_riming_snow[i, k] * ratio
            qc_bergeron_snow[i, k] = qc_bergeron_snow[i, k] * ratio
            qc_bergeron_ice[i, k] = qc_bergeron_ice[i, k] * ratio
            qc_limiter_ratio[i, k] = ratio
        else
            qc_limiter_ratio[i, k] = 1.0
        end
        if qc[i, k] >= smallest_mixing_ratio
            qi_deposition[i, k] = qi_deposition[i, k] * (1.0 - qc_limiter_ratio[i, k])
            qs_deposition[i, k] = qs_deposition[i, k] * (1.0 - qc_limiter_ratio[i, k])
        end
        qi_nucleation[i, k] = max(0.0, qi_nucleation[i, k])
        qi_deposition[i, k] = max(0.0, qi_deposition[i, k])
        qs_deposition[i, k] = max(0.0, qs_deposition[i, k])
        deposition_total = qi_deposition[i, k] + qi_nucleation[i, k] + qs_deposition[i, k]
        if deposition_total > 1.0e-20
            deposition_max = (qv[i, k] - qsat_ice[i, k]) /
                             (1.0 +
                              Ls_squared * qsat_ice[i, k] /
                              (cp * Rv * (T[i, k] * T[i, k]))) *
                             inv_Δt
            deposition_max = max(deposition_max, 0.0)
            if deposition_total > deposition_max
                qi_nucleation[i, k] = deposition_max * qi_nucleation[i, k] /
                                      deposition_total
                qi_deposition[i, k] = deposition_max * qi_deposition[i, k] /
                                      deposition_total
                qs_deposition[i, k] = deposition_max * qs_deposition[i, k] /
                                      deposition_total
            end
        end
    end

    # 💧🔢 conservation of nc; small raindrops freeze to cloud ice
    for k in 1:levels, i in 1:columns
        sink_over_step = (nc_autoconversion[i, k] + nc_accretion[i, k] +
                          nc_immersion_freezing[i, k] + nc_contact_freezing[i, k] +
                          nc_riming_snow[i, k] -
                          nc_evaporation[i, k] + nc_riming_graupel[i, k]) *
                         liquid_fraction[i, k] * Δt
        if sink_over_step > nc[i, k]
            ratio = nc[i, k] * inv_Δt /
                    ((nc_autoconversion[i, k] + nc_accretion[i, k] +
                      nc_immersion_freezing[i, k] + nc_contact_freezing[i, k] +
                      nc_riming_snow[i, k] - nc_evaporation[i, k] +
                      nc_riming_graupel[i, k]) * liquid_fraction[i, k]) * one_minus_small
            nc_riming_graupel[i, k] = nc_riming_graupel[i, k] * ratio
            nc_autoconversion[i, k] = nc_autoconversion[i, k] * ratio
            nc_accretion[i, k] = nc_accretion[i, k] * ratio
            nc_immersion_freezing[i, k] = nc_immersion_freezing[i, k] * ratio
            nc_contact_freezing[i, k] = nc_contact_freezing[i, k] * ratio
            nc_riming_snow[i, k] = nc_riming_snow[i, k] * ratio
            nc_evaporation[i, k] = nc_evaporation[i, k] * ratio
        end
        qr_freezing_ice[i, k] = 0.0
        nr_freezing_ice[i, k] = 0.0
        if settings.prognostic_cloud_ice && λr[i, k] > smallest_mixing_ratio &&
           1.0 / λr[i, k] < settings.ice_to_snow_size
            qr_freezing_ice[i, k] = qr_freezing_snow[i, k]
            nr_freezing_ice[i, k] = nr_freezing_snow[i, k]
            qr_freezing_snow[i, k] = 0.0
            nr_freezing_snow[i, k] = 0.0
        end
    end

    # 🌧️ conservation of rain mass and number
    for k in 1:levels, i in 1:columns
        sink_over_step = ((-qr_evaporation[i, k] + qr_accretion_snow[i, k] +
                           qr_freezing_snow[i, k] + qr_freezing_ice[i, k] +
                           qi_splintering_rain[i, k] +
                           qr_accretion_graupel[i, k] + qg_rain_snow[i, k]) *
                          precip_fraction[i, k] -
                          (qc_accretion[i, k] + qc_autoconversion[i, k]) *
                          liquid_fraction[i, k]) * Δt
        sinks = -qr_evaporation[i, k] + qr_accretion_snow[i, k] + qr_freezing_snow[i, k] +
                qr_freezing_ice[i, k] + qi_splintering_rain[i, k] +
                qr_accretion_graupel[i, k] + qg_rain_snow[i, k]
        if sink_over_step > qr[i, k] && sinks >= smallest_mixing_ratio
            ratio = (qr[i, k] * inv_Δt +
                     (qc_accretion[i, k] + qc_autoconversion[i, k]) *
                     liquid_fraction[i, k]) /
                    precip_fraction[i, k] / sinks * one_minus_small
            qi_splintering_rain[i, k] = qi_splintering_rain[i, k] * ratio
            qr_accretion_graupel[i, k] = qr_accretion_graupel[i, k] * ratio
            qg_rain_snow[i, k] = qg_rain_snow[i, k] * ratio
            qr_evaporation[i, k] = qr_evaporation[i, k] * ratio
            qr_accretion_snow[i, k] = qr_accretion_snow[i, k] * ratio
            qr_freezing_snow[i, k] = qr_freezing_snow[i, k] * ratio
            qr_freezing_ice[i, k] = qr_freezing_ice[i, k] * ratio
        end
        if qr_evaporation[i, k] < 0.0
            nr_evaporation[i, k] = qr_evaporation[i, k] * nr[i, k] / qr[i, k]
        else
            nr_evaporation[i, k] = 0.0
        end
        nr_sinks = -nr_evaporation[i, k] + nr_accretion_snow[i, k] +
                   nr_freezing_snow[i, k] + nr_freezing_ice[i, k] -
                   nr_self_collection[i, k] +
                   nr_accretion_graupel[i, k] + ng_rain_snow[i, k]
        sink_over_step = (nr_sinks * precip_fraction[i, k] -
                          nr_autoconversion[i, k] * liquid_fraction[i, k]) * Δt
        nr_sources = nr[i, k] * inv_Δt + nr_autoconversion[i, k] * liquid_fraction[i, k]
        if sink_over_step > nr[i, k] && nr_sinks > 0.0 && nr_sources > 0.0 &&
           precip_fraction[i, k] > 0.0
            ratio = (nr[i, k] * inv_Δt + nr_autoconversion[i, k] * liquid_fraction[i, k]) /
                    precip_fraction[i, k] /
                    nr_sinks * one_minus_small
            nr_accretion_graupel[i, k] = nr_accretion_graupel[i, k] * ratio
            ng_rain_snow[i, k] = ng_rain_snow[i, k] * ratio
            nr_self_collection[i, k] = nr_self_collection[i, k] * ratio
            nr_accretion_snow[i, k] = nr_accretion_snow[i, k] * ratio
            nr_freezing_snow[i, k] = nr_freezing_snow[i, k] * ratio
            nr_evaporation[i, k] = nr_evaporation[i, k] * ratio
            nr_freezing_ice[i, k] = nr_freezing_ice[i, k] * ratio
        end
    end

    # 🧊 conservation of cloud ice mass and number
    if settings.prognostic_cloud_ice
        for k in 1:levels, i in 1:columns
            sink_over_step = ((-qc_immersion_freezing[i, k] - qc_contact_freezing[i, k] -
                               qc_deposition_nucleation[i, k] - qi_splintering_snow[i, k] -
                               qi_splintering_graupel[i, k]) * liquid_fraction[i, k] +
                              (qi_autoconversion[i, k] + qi_accretion_snow[i, k]) *
                              ice_fraction[i, k] +
                              (-qi_splintering_rain[i, k] - qr_freezing_ice[i, k]) *
                              precip_fraction[i, k] - qi_sublimation[i, k] -
                              qi_deposition[i, k] - qc_bergeron_ice[i, k] -
                              qi_nucleation[i, k]) * Δt
            if sink_over_step > qi[i, k]
                ratio = (qi[i, k] * inv_Δt + qi_deposition[i, k] + qc_bergeron_ice[i, k] +
                         qi_nucleation[i, k] +
                         (qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                          qc_deposition_nucleation[i, k] + qi_splintering_snow[i, k] +
                          qi_splintering_graupel[i, k]) * liquid_fraction[i, k] +
                         (qi_splintering_rain[i, k] + qr_freezing_ice[i, k]) *
                         precip_fraction[i, k]) /
                        ((qi_autoconversion[i, k] + qi_accretion_snow[i, k]) *
                         ice_fraction[i, k] - qi_sublimation[i, k]) * one_minus_small
                qi_autoconversion[i, k] = qi_autoconversion[i, k] * ratio
                qi_accretion_snow[i, k] = qi_accretion_snow[i, k] * ratio
                qi_sublimation[i, k] = qi_sublimation[i, k] * ratio
            end
            nc_immersion_counted = settings.host_heterogeneous_freezing ?
                                   nc_immersion_freezing[i, k] : 0.0
            sink_over_step = ((-nc_contact_freezing[i, k] - nc_immersion_counted -
                               nc_deposition_nucleation[i, k] - ni_splintering_snow[i, k] -
                               ni_splintering_graupel[i, k]) *
                              liquid_fraction[i, k] +
                              (ni_autoconversion[i, k] + ni_accretion_snow[i, k] -
                               ni_sublimation[i, k]) * ice_fraction[i, k] +
                              (-ni_splintering_rain[i, k] - nr_freezing_ice[i, k]) *
                              precip_fraction[i, k] - ni_nucleation[i, k]) *
                             Δt
            if sink_over_step > ni[i, k]
                ratio = (ni[i, k] * inv_Δt + ni_nucleation[i, k] +
                         (nc_contact_freezing[i, k] + nc_immersion_counted +
                          nc_deposition_nucleation[i, k] + ni_splintering_snow[i, k] +
                          ni_splintering_graupel[i, k]) * liquid_fraction[i, k] +
                         (nr_freezing_ice[i, k] + ni_splintering_rain[i, k]) *
                         precip_fraction[i, k]) /
                        ((ni_autoconversion[i, k] + ni_accretion_snow[i, k] -
                          ni_sublimation[i, k]) * ice_fraction[i, k]) * one_minus_small
                ni_autoconversion[i, k] = ni_autoconversion[i, k] * ratio
                ni_accretion_snow[i, k] = ni_accretion_snow[i, k] * ratio
                ni_sublimation[i, k] = ni_sublimation[i, k] * ratio
            end
        end
    end

    # ❄️ conservation of snow mass and number
    for k in 1:levels, i in 1:columns
        if graupel_on
            sink_over_step = (-(qs_sublimation[i, k] + qr_accretion_snow[i, k] -
                                qs_collection_rain[i, k]) * precip_fraction[i, k] -
                              (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) *
                              ice_fraction[i, k] -
                              (qc_bergeron_snow[i, k] + qc_riming_snow[i, k]) *
                              liquid_fraction[i, k] - qs_deposition[i, k]) * Δt
        else
            sink_over_step = (-(qs_sublimation[i, k] + qr_accretion_snow[i, k] +
                                qr_freezing_snow[i, k]) * precip_fraction[i, k] -
                              (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) *
                              ice_fraction[i, k] -
                              (qc_bergeron_snow[i, k] + qc_riming_snow[i, k]) *
                              liquid_fraction[i, k] - qs_deposition[i, k]) * Δt
        end
        if sink_over_step > qs[i, k] &&
           (qs_collection_rain[i, k] - qs_sublimation[i, k]) >= smallest_mixing_ratio
            if graupel_on
                ratio = (qs[i, k] * inv_Δt +
                         (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) *
                         ice_fraction[i, k] +
                         (qc_bergeron_snow[i, k] + qc_riming_snow[i, k]) *
                         liquid_fraction[i, k] + qs_deposition[i, k] +
                         qr_accretion_snow[i, k] * precip_fraction[i, k]) /
                        precip_fraction[i, k] /
                        (qs_collection_rain[i, k] - qs_sublimation[i, k]) * one_minus_small
                qs_collection_rain[i, k] = qs_collection_rain[i, k] * ratio
            else
                ratio = (qs[i, k] * inv_Δt +
                         (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) *
                         ice_fraction[i, k] +
                         (qc_bergeron_snow[i, k] + qc_riming_snow[i, k]) *
                         liquid_fraction[i, k] + qs_deposition[i, k] +
                         (qr_accretion_snow[i, k] + qr_freezing_snow[i, k]) *
                         precip_fraction[i, k]) /
                        precip_fraction[i, k] / (-qs_sublimation[i, k]) * one_minus_small
            end
            qs_sublimation[i, k] = qs_sublimation[i, k] * ratio
        end
        ns_sublimation[i, k] = 0.0
        if graupel_on
            sink_over_step = ((-ns_aggregation[i, k] - ns_sublimation[i, k] +
                               ng_rain_snow[i, k]) * precip_fraction[i, k] -
                              ni_autoconversion[i, k] * ice_fraction[i, k] +
                              ng_rimed_snow[i, k] * liquid_fraction[i, k]) * Δt
        else
            sink_over_step = ((-ns_aggregation[i, k] - ns_sublimation[i, k] -
                               nr_freezing_snow[i, k]) * precip_fraction[i, k] -
                              ni_autoconversion[i, k] * ice_fraction[i, k]) * Δt
        end
        if sink_over_step > ns[i, k]
            if graupel_on
                ratio = (ns[i, k] * inv_Δt + ni_autoconversion[i, k] * ice_fraction[i, k]) /
                        precip_fraction[i, k] /
                        (-ns_sublimation[i, k] - ns_aggregation[i, k] + ng_rain_snow[i, k] +
                         liquid_fraction[i, k] / precip_fraction[i, k] *
                         ng_rimed_snow[i, k]) *
                        one_minus_small
                ng_rimed_snow[i, k] = ng_rimed_snow[i, k] * ratio
                ng_rain_snow[i, k] = ng_rain_snow[i, k] * ratio
            else
                ratio = (ns[i, k] * inv_Δt +
                         nr_freezing_snow[i, k] * precip_fraction[i, k] +
                         ni_autoconversion[i, k] * ice_fraction[i, k]) /
                        precip_fraction[i, k] /
                        (-ns_sublimation[i, k] - ns_aggregation[i, k]) * one_minus_small
            end
            ns_sublimation[i, k] = ns_sublimation[i, k] * ratio
            ns_aggregation[i, k] = ns_aggregation[i, k] * ratio
        end
    end

    # 🟤 conservation of graupel mass (graupel number has no sinks)
    if graupel_on
        for k in 1:levels, i in 1:columns
            sink_over_step = ((-qr_accretion_graupel[i, k] - qg_rain_snow[i, k] -
                               qg_deposition[i, k] - qs_collection_rain[i, k] -
                               qr_freezing_snow[i, k]) *
                              precip_fraction[i, k] +
                              (-qc_riming_graupel[i, k] - qg_rimed_snow[i, k]) *
                              liquid_fraction[i, k]) * Δt
            if sink_over_step > qg[i, k]
                ratio = (qg[i, k] * inv_Δt +
                         (qr_accretion_graupel[i, k] + qg_rain_snow[i, k] +
                          qs_collection_rain[i, k] + qr_freezing_snow[i, k]) *
                         precip_fraction[i, k] +
                         (qc_riming_graupel[i, k] + qg_rimed_snow[i, k]) *
                         liquid_fraction[i, k]) /
                        ((-qg_deposition[i, k]) * precip_fraction[i, k]) * one_minus_small
                qg_deposition[i, k] = qg_deposition[i, k] * ratio
            end
        end
    end

    # 🌫️ limit evaporation and sublimation so that they do not supersaturate the air
    for k in 1:levels, i in 1:columns
        qv_estimate[i, k] = qv[i, k] -
                            (qi_sublimation[i, k] + qi_deposition[i, k] +
                             qi_nucleation[i, k] +
                             qs_deposition[i, k] +
                             (qr_evaporation[i, k] + qs_sublimation[i, k] +
                              qg_deposition[i, k]) * precip_fraction[i, k]) *
                            Δt
        T_estimate[i, k] = T[i, k] +
                           ((qr_evaporation[i, k] * precip_fraction[i, k]) * Lv +
                            ((qs_sublimation[i, k] + qg_deposition[i, k]) *
                             precip_fraction[i, k] +
                             qi_deposition[i, k] + qs_deposition[i, k] +
                             qi_sublimation[i, k] +
                             qi_nucleation[i, k]) * Ls) * Δt / cp
    end
    for k in 1:levels, i in 1:columns
        qsat_estimate[i, k] = saturation_over_water(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    function evaporating(i, k)
        (qr_evaporation[i, k] + qs_sublimation[i, k] + qg_deposition[i, k]) *
        precip_fraction[i, k] + qi_sublimation[i, k] <
        -1.0e-20 && qv_estimate[i, k] > qsat_estimate[i, k]
    end
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            total = (qr_evaporation[i, k] + qs_sublimation[i, k] + qg_deposition[i, k]) *
                    precip_fraction[i, k] +
                    qi_sublimation[i, k]
            rain_share[i, k] = qr_evaporation[i, k] * precip_fraction[i, k] / total
            snow_share[i, k] = qs_sublimation[i, k] * precip_fraction[i, k] / total
            graupel_share[i, k] = qg_deposition[i, k] * precip_fraction[i, k] / total
            T_estimate[i, k] = T[i, k] +
                               ((qi_deposition[i, k] + qs_deposition[i, k] +
                                 qi_nucleation[i, k]) * Ls) *
                               Δt / cp
            qv_after_deposition[i, k] = qv[i, k] -
                                        (qi_deposition[i, k] + qs_deposition[i, k] +
                                         qi_nucleation[i, k]) * Δt
        end
    end
    for k in 1:levels, i in 1:columns
        qsat_adjusted[i, k] = saturation_over_water(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            vapor_excess = (qv_after_deposition[i, k] - qsat_adjusted[i, k]) /
                           (1.0 +
                            Lv_squared * qsat_adjusted[i, k] /
                            (cp * Rv * (T_estimate[i, k] * T_estimate[i, k])))
            vapor_excess = min(vapor_excess, 0.0)
            qr_evaporation[i, k] = vapor_excess * rain_share[i, k] * inv_Δt /
                                   precip_fraction[i, k]
        end
    end
    for k in 1:levels, i in 1:columns
        qsat_adjusted[i, k] = saturation_over_ice(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        if evaporating(i, k)
            vapor_excess = (qv_after_deposition[i, k] - qsat_adjusted[i, k]) /
                           (1.0 +
                            Ls_squared * qsat_adjusted[i, k] /
                            (cp * Rv * (T_estimate[i, k] * T_estimate[i, k])))
            vapor_excess = min(vapor_excess, 0.0)
            qs_sublimation[i, k] = vapor_excess * snow_share[i, k] * inv_Δt /
                                   precip_fraction[i, k]
            qg_deposition[i, k] = vapor_excess * graupel_share[i, k] * inv_Δt /
                                  precip_fraction[i, k]
            # 🧊 rain_share (dum1A) is reused for the cloud ice share
            rain_share[i, k] = (1.0 - rain_share[i, k] - snow_share[i, k] -
                                graupel_share[i, k])
            qi_sublimation[i, k] = vapor_excess * rain_share[i, k] * inv_Δt
        end
    end
    return nothing
end
