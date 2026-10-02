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
    (; Lv, Ls, Lf) = host_constants(parameters)
    (; state, incloud, rates, process_rates, tendencies, diagnostics) = work
    qi_net_deposition = process_rates.cloud_ice.net_deposition.mass  # cmeitot
    qs_sublimation_mean = process_rates.snow.sublimation.mass  # evapsnow
    qs_production = process_rates.snow.production.mass  # prodsnow
    qc_accretion_mean = process_rates.cloud_liquid.accretion_by_rain.mass  # pratot
    qc_autoconversion_mean = process_rates.cloud_liquid.autoconversion_to_rain.mass  # prctot
    qc_immersion_freezing_mean = process_rates.cloud_liquid.immersion_freezing.mass  # mnuccctot
    qc_deposition_nucleation_mean = process_rates.cloud_liquid.deposition_nucleation.mass  # mnudeptot
    qc_contact_freezing_mean = process_rates.cloud_liquid.contact_freezing.mass  # mnuccttot
    qi_splintering_snow_mean = process_rates.cloud_ice.splintering_by_snow.mass  # msacwitot
    qc_riming_snow_mean = process_rates.cloud_liquid.riming_by_snow.mass  # psacwstot
    qc_bergeron_snow_mean = process_rates.cloud_liquid.bergeron_to_snow.mass  # bergstot
    qs_deposition_mean = process_rates.snow.deposition.mass  # vapdepstot
    qc_bergeron_ice_mean = process_rates.cloud_liquid.bergeron_to_cloud_ice.mass  # bergtot
    qi_autoconversion_mean = process_rates.cloud_ice.autoconversion_to_snow.mass  # prcitot
    qi_accretion_snow_mean = process_rates.cloud_ice.accretion_by_snow.mass  # praitot
    qi_nucleation_mean = process_rates.cloud_ice.nucleation_from_vapor.mass  # mnuccdtot
    qr_accretion_snow_mean = process_rates.rain.accretion_by_snow.mass  # pracstot
    qr_freezing_snow_mean = process_rates.rain.freezing_to_snow.mass  # mnuccrtot
    qr_freezing_ice_mean = process_rates.rain.freezing_to_cloud_ice.mass  # mnuccritot
    qs_collection_rain_mean = process_rates.snow.collection_by_rain.mass  # psacrtot
    qr_accretion_graupel_mean = process_rates.rain.accretion_by_graupel.mass  # pracgtot
    qc_riming_graupel_mean = process_rates.cloud_liquid.riming_by_graupel.mass  # psacwgtot
    qg_rimed_snow_mean = process_rates.graupel.from_snow_riming.mass  # pgsacwtot
    qg_rain_snow_mean = process_rates.graupel.rain_snow_collection.mass  # pgracstot
    qg_deposition_mean = process_rates.graupel.deposition.mass  # prdgtot
    qi_splintering_graupel_mean = process_rates.cloud_ice.splintering_by_graupel.mass  # qmultgtot
    qi_splintering_rain_mean = process_rates.cloud_ice.splintering_rain_graupel.mass  # qmultrgtot
    nr_accretion_graupel_mean = process_rates.rain.accretion_by_graupel.number  # npracgtot
    ng_rimed_snow_mean = process_rates.graupel.from_snow_riming.number  # nscngtot
    ng_rain_snow_mean = process_rates.graupel.rain_snow_collection.number  # ngracstot
    ni_splintering_graupel_mean = process_rates.cloud_ice.splintering_by_graupel.number  # nmultgtot
    ni_splintering_rain_mean = process_rates.cloud_ice.splintering_rain_graupel.number  # nmultrgtot
    nc_riming_graupel_mean = process_rates.cloud_liquid.riming_by_graupel.number  # npsacwgtot
    nc_immersion_freezing_mean = process_rates.cloud_liquid.immersion_freezing.number  # nnuccctot
    nc_contact_freezing_mean = process_rates.cloud_liquid.contact_freezing.number  # nnuccttot
    ni_nucleation_mean = process_rates.cloud_ice.nucleation_from_vapor.number  # nnuccdtot
    nc_deposition_nucleation_mean = process_rates.cloud_liquid.deposition_nucleation.number  # nnudeptot
    nr_freezing_snow_mean = process_rates.rain.freezing_to_snow.number  # nnuccrtot
    nr_freezing_ice_mean = process_rates.rain.freezing_to_cloud_ice.number  # nnuccritot
    ni_splintering_snow_mean = process_rates.cloud_ice.splintering_by_snow.number  # nsacwitot
    nc_accretion_mean = process_rates.cloud_liquid.accretion_by_rain.number  # npratot
    nc_riming_snow_mean = process_rates.cloud_liquid.riming_by_snow.number  # npsacwstot
    ni_accretion_snow_mean = process_rates.cloud_ice.accretion_by_snow.number  # npraitot
    nr_accretion_snow_mean = process_rates.rain.accretion_by_snow.number  # npracstot
    nr_autoconversion_mean = process_rates.rain.autoconversion_from_cloud_liquid.number  # nprctot
    nr_self_collection_mean = process_rates.rain.self_collection.number  # nraggtot
    ni_autoconversion_mean = process_rates.cloud_ice.autoconversion_to_snow.number  # nprcitot
    ns_melting_mean = process_rates.snow.melting.number  # nmeltstot
    ng_melting_mean = process_rates.graupel.melting.number  # nmeltgtot
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    qr = state.rain.q
    nr = state.rain.N
    qs = state.snow.q
    ns = state.snow.N
    qg = state.graupel.q
    ng = state.graupel.N
    ρ = state.air.density  # rho
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    qv_tendency = tendencies.vapor  # qvlat
    heating = tendencies.heating  # tlat
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
    ns_instant_melting = rates.snow.instant_melting.number  # ninstsm
    ng_instant_melting = rates.graupel.instant_melting.number  # ninstgm
    net_condensation = diagnostics.vapor.net_condensation  # cmeout
    precip_evaporation = diagnostics.precipitation.evaporation  # nevapr
    rain_evaporation = diagnostics.rain.evaporation  # prer_evap
    precip_production = diagnostics.precipitation.production  # prain
    qc_sink_rate = diagnostics.cloud_liquid.first_order_sink_rate  # qcsinksum_rate1ord
    qr_out = diagnostics.rain.q  # qrout
    nr_per_volume = diagnostics.rain.N_per_volume  # nrout
    qs_out = diagnostics.snow.q  # qsout
    ns_per_volume = diagnostics.snow.N_per_volume  # nsout
    qg_out = diagnostics.graupel.q  # qgout
    ng_per_volume = diagnostics.graupel.N_per_volume  # ngout
    qr_precip = incloud.rain.q  # qric
    nr_precip = incloud.rain.N  # nric
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    λc = incloud.cloud_liquid.λ  # lamc
    μc = incloud.cloud_liquid.μ  # pgam
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    radius_rain_and_cloud = diagnostics.effective_radius.rain_and_cloud  # rercld
    Δt = work.inputs.time_step
    graupel_on = settings.do_hail || settings.do_graupel
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        qv_tendency[i, k] = qv_tendency[i, k] - (qr_evaporation[i, k] + qs_sublimation[i, k]) * precip_fraction[i, k] -
                      qi_deposition[i, k] - qs_deposition[i, k] - qi_sublimation[i, k] - qi_nucleation[i, k] -
                      qc_deposition_nucleation[i, k] * liquid_fraction[i, k] - qg_deposition[i, k] * precip_fraction[i, k]
        heating[i, k] = heating[i, k] +
                     ((qr_evaporation[i, k] * precip_fraction[i, k]) * Lv +
                      ((qs_sublimation[i, k] + qg_deposition[i, k]) * precip_fraction[i, k] + qi_deposition[i, k] +
                       qs_deposition[i, k] + qi_sublimation[i, k] + qi_nucleation[i, k] +
                       qc_deposition_nucleation[i, k] * liquid_fraction[i, k]) * Ls +
                      ((qc_bergeron_snow[i, k] + qc_riming_snow[i, k] + qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                        qi_splintering_snow[i, k] + qc_riming_graupel[i, k] + qi_splintering_graupel[i, k] + qg_rimed_snow[i, k]) *
                       liquid_fraction[i, k] +
                       (qr_freezing_snow[i, k] + qr_accretion_snow[i, k] + qr_freezing_ice[i, k] + qr_accretion_graupel[i, k] +
                        qg_rain_snow[i, k] + qi_splintering_rain[i, k]) * precip_fraction[i, k] +
                       qc_bergeron_ice[i, k]) * Lf)
        qc_tendency[i, k] = qc_tendency[i, k] +
                       (-qc_accretion[i, k] - qc_autoconversion[i, k] - qc_immersion_freezing[i, k] - qc_contact_freezing[i, k] -
                        qi_splintering_snow[i, k] - qc_riming_snow[i, k] - qc_bergeron_snow[i, k] - qi_splintering_graupel[i, k] -
                        qc_riming_graupel[i, k] - qg_rimed_snow[i, k]) * liquid_fraction[i, k] - qc_bergeron_ice[i, k]
        if settings.prognostic_cloud_ice
            qi_tendency[i, k] = qi_tendency[i, k] +
                           (qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] + qc_deposition_nucleation[i, k] + qi_splintering_snow[i, k] +
                            qi_splintering_graupel[i, k]) * liquid_fraction[i, k] +
                           (-qi_autoconversion[i, k] - qi_accretion_snow[i, k]) * ice_fraction[i, k] + qi_deposition[i, k] +
                           qc_bergeron_ice[i, k] + qi_sublimation[i, k] + qi_nucleation[i, k] +
                           (qr_freezing_ice[i, k] + qi_splintering_rain[i, k]) * precip_fraction[i, k]
        end
        qr_tendency[i, k] = qr_tendency[i, k] + (qc_accretion[i, k] + qc_autoconversion[i, k]) * liquid_fraction[i, k] +
                       (qr_evaporation[i, k] - qr_accretion_snow[i, k] - qr_freezing_snow[i, k] - qr_freezing_ice[i, k] -
                        qi_splintering_rain[i, k] - qr_accretion_graupel[i, k] - qg_rain_snow[i, k]) * precip_fraction[i, k]
        if graupel_on
            qg_tendency[i, k] = qg_tendency[i, k] +
                           (qr_accretion_graupel[i, k] + qg_rain_snow[i, k] + qg_deposition[i, k] + qs_collection_rain[i, k] +
                            qr_freezing_snow[i, k]) * precip_fraction[i, k] +
                           (qc_riming_graupel[i, k] + qg_rimed_snow[i, k]) * liquid_fraction[i, k]
            qs_tendency[i, k] = qs_tendency[i, k] + (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) * ice_fraction[i, k] +
                           (qc_riming_snow[i, k] + qc_bergeron_snow[i, k]) * liquid_fraction[i, k] +
                           (qs_sublimation[i, k] + qr_accretion_snow[i, k] - qs_collection_rain[i, k]) * precip_fraction[i, k] +
                           qs_deposition[i, k]
        else
            qs_tendency[i, k] = qs_tendency[i, k] + (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) * ice_fraction[i, k] +
                           (qc_riming_snow[i, k] + qc_bergeron_snow[i, k]) * liquid_fraction[i, k] +
                           (qs_sublimation[i, k] + qr_accretion_snow[i, k] + qr_freezing_snow[i, k]) * precip_fraction[i, k] +
                           qs_deposition[i, k]
        end
    end

    for k in 1:levels, i in 1:columns
        net_condensation[i, k] = qi_deposition[i, k] + qi_sublimation[i, k] + qi_nucleation[i, k] + qs_deposition[i, k]
        qi_net_deposition[i, k] = qi_deposition[i, k] + qi_sublimation[i, k] + qi_nucleation[i, k] + qs_deposition[i, k]
        qs_sublimation_mean[i, k] = (-qs_sublimation[i, k] - qg_deposition[i, k]) * precip_fraction[i, k]
        precip_evaporation[i, k] = -qr_evaporation[i, k] * precip_fraction[i, k]
        rain_evaporation[i, k] = -qr_evaporation[i, k] * precip_fraction[i, k]
        precip_production[i, k] = (qc_accretion[i, k] + qc_autoconversion[i, k]) * liquid_fraction[i, k] +
                      (-qr_accretion_snow[i, k] - qr_freezing_snow[i, k] - qr_freezing_ice[i, k]) * precip_fraction[i, k]
        if graupel_on
            qs_production[i, k] = (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) * ice_fraction[i, k] +
                               (qc_riming_snow[i, k] + qc_bergeron_snow[i, k]) * liquid_fraction[i, k] +
                               (qr_accretion_snow[i, k]) * precip_fraction[i, k] + qs_deposition[i, k]
        else
            qs_production[i, k] = (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) * ice_fraction[i, k] +
                               (qc_riming_snow[i, k] + qc_bergeron_snow[i, k]) * liquid_fraction[i, k] +
                               (qr_accretion_snow[i, k] + qr_freezing_snow[i, k]) * precip_fraction[i, k] +
                               qs_deposition[i, k]
        end
        qc_sink_rate[i, k] = (qc_accretion[i, k] + qc_autoconversion[i, k] + qc_riming_snow[i, k] + qc_riming_graupel[i, k] +
                                    qg_rimed_snow[i, k]) * liquid_fraction[i, k]
        qc_sink_rate[i, k] = qc_sink_rate[i, k] / max(qc[i, k], 1.0e-30)
    end

    for k in 1:levels, i in 1:columns
        qc_accretion_mean[i, k] = qc_accretion[i, k] * liquid_fraction[i, k]
        qc_autoconversion_mean[i, k] = qc_autoconversion[i, k] * liquid_fraction[i, k]
        qc_immersion_freezing_mean[i, k] = qc_immersion_freezing[i, k] * liquid_fraction[i, k]
        qc_deposition_nucleation_mean[i, k] = qc_deposition_nucleation[i, k] * liquid_fraction[i, k]
        qc_contact_freezing_mean[i, k] = qc_contact_freezing[i, k] * liquid_fraction[i, k]
        qi_splintering_snow_mean[i, k] = qi_splintering_snow[i, k] * liquid_fraction[i, k]
        qc_riming_snow_mean[i, k] = qc_riming_snow[i, k] * liquid_fraction[i, k]
        qc_bergeron_snow_mean[i, k] = qc_bergeron_snow[i, k] * liquid_fraction[i, k]
        qs_deposition_mean[i, k] = qs_deposition[i, k]
        qc_bergeron_ice_mean[i, k] = qc_bergeron_ice[i, k]
        qi_autoconversion_mean[i, k] = qi_autoconversion[i, k] * ice_fraction[i, k]
        qi_accretion_snow_mean[i, k] = qi_accretion_snow[i, k] * ice_fraction[i, k]
        qi_nucleation_mean[i, k] = qi_nucleation[i, k] * ice_fraction[i, k]
        qr_accretion_snow_mean[i, k] = qr_accretion_snow[i, k] * precip_fraction[i, k]
        qr_freezing_snow_mean[i, k] = qr_freezing_snow[i, k] * precip_fraction[i, k]
        qr_freezing_ice_mean[i, k] = qr_freezing_ice[i, k] * precip_fraction[i, k]
        qs_collection_rain_mean[i, k] = qs_collection_rain[i, k] * precip_fraction[i, k]
        qr_accretion_graupel_mean[i, k] = qr_accretion_graupel[i, k] * precip_fraction[i, k]
        qc_riming_graupel_mean[i, k] = qc_riming_graupel[i, k] * liquid_fraction[i, k]
        qg_rimed_snow_mean[i, k] = qg_rimed_snow[i, k] * liquid_fraction[i, k]
        qg_rain_snow_mean[i, k] = qg_rain_snow[i, k] * precip_fraction[i, k]
        qg_deposition_mean[i, k] = qg_deposition[i, k] * precip_fraction[i, k]
        qi_splintering_graupel_mean[i, k] = qi_splintering_graupel[i, k] * liquid_fraction[i, k]
        qi_splintering_rain_mean[i, k] = qi_splintering_rain[i, k] * precip_fraction[i, k]
        nr_accretion_graupel_mean[i, k] = nr_accretion_graupel[i, k] * precip_fraction[i, k]
        ng_rimed_snow_mean[i, k] = ng_rimed_snow[i, k] * liquid_fraction[i, k]
        ng_rain_snow_mean[i, k] = ng_rain_snow[i, k] * precip_fraction[i, k]
        ni_splintering_graupel_mean[i, k] = ni_splintering_graupel[i, k] * liquid_fraction[i, k]
        ni_splintering_rain_mean[i, k] = ni_splintering_rain[i, k] * precip_fraction[i, k]
        nc_riming_graupel_mean[i, k] = nc_riming_graupel[i, k] * liquid_fraction[i, k]
        nc_immersion_freezing_mean[i, k] = nc_immersion_freezing[i, k] * liquid_fraction[i, k]
        nc_contact_freezing_mean[i, k] = nc_contact_freezing[i, k] * liquid_fraction[i, k]
        ni_nucleation_mean[i, k] = ni_nucleation[i, k] * ice_fraction[i, k]
        nc_deposition_nucleation_mean[i, k] = nc_deposition_nucleation[i, k] * liquid_fraction[i, k]
        nr_freezing_snow_mean[i, k] = nr_freezing_snow[i, k] * precip_fraction[i, k]
        nr_freezing_ice_mean[i, k] = nr_freezing_ice[i, k] * precip_fraction[i, k]
        ni_splintering_snow_mean[i, k] = ni_splintering_snow[i, k] * liquid_fraction[i, k]
        nc_accretion_mean[i, k] = nc_accretion[i, k] * liquid_fraction[i, k]
        nc_riming_snow_mean[i, k] = nc_riming_snow[i, k] * liquid_fraction[i, k]
        ni_accretion_snow_mean[i, k] = ni_accretion_snow[i, k] * ice_fraction[i, k]
        nr_accretion_snow_mean[i, k] = nr_accretion_snow[i, k] * precip_fraction[i, k]
        nr_autoconversion_mean[i, k] = nr_autoconversion[i, k] * liquid_fraction[i, k]
        nr_self_collection_mean[i, k] = nr_self_collection[i, k] * precip_fraction[i, k]
        ni_autoconversion_mean[i, k] = ni_autoconversion[i, k] * ice_fraction[i, k]
        ns_melting_mean[i, k] = ns_instant_melting[i, k] / Δt
        ng_melting_mean[i, k] = ng_instant_melting[i, k] / Δt
    end

    for k in 1:levels, i in 1:columns
        nc_tendency[i, k] = nc_tendency[i, k] +
                       (-nc_immersion_freezing[i, k] - nc_contact_freezing[i, k] - nc_riming_snow[i, k] + nc_evaporation[i, k] -
                        nc_accretion[i, k] - nc_autoconversion[i, k] - nc_riming_graupel[i, k]) * liquid_fraction[i, k]
        if settings.prognostic_cloud_ice
            tmpfrz = settings.host_heterogeneous_freezing ? nc_immersion_freezing[i, k] : 0.0
            ni_tendency[i, k] = ni_tendency[i, k] + ni_nucleation[i, k] +
                           (nc_contact_freezing[i, k] + tmpfrz + nc_deposition_nucleation[i, k] + ni_splintering_snow[i, k] +
                            ni_splintering_graupel[i, k]) * liquid_fraction[i, k] +
                           (ni_sublimation[i, k] - ni_autoconversion[i, k] - ni_accretion_snow[i, k]) * ice_fraction[i, k] +
                           (nr_freezing_ice[i, k] + ni_splintering_rain[i, k]) * precip_fraction[i, k]
        end
        if graupel_on
            ns_tendency[i, k] = ns_tendency[i, k] +
                           (ns_sublimation[i, k] + ns_aggregation[i, k] - ng_rain_snow[i, k]) * precip_fraction[i, k] +
                           ni_autoconversion[i, k] * ice_fraction[i, k] - ng_rimed_snow[i, k] * liquid_fraction[i, k]
            ng_tendency[i, k] = ng_tendency[i, k] + ng_rimed_snow[i, k] * liquid_fraction[i, k] +
                           (ng_rain_snow[i, k] + nr_freezing_snow[i, k]) * precip_fraction[i, k]
        else
            ns_tendency[i, k] = ns_tendency[i, k] +
                           (ns_sublimation[i, k] + ns_aggregation[i, k] + nr_freezing_snow[i, k]) * precip_fraction[i, k] +
                           ni_autoconversion[i, k] * ice_fraction[i, k]
        end
        nr_tendency[i, k] = nr_tendency[i, k] + nr_autoconversion[i, k] * liquid_fraction[i, k] +
                       (nr_evaporation[i, k] - nr_accretion_snow[i, k] - nr_freezing_snow[i, k] - nr_freezing_ice[i, k] +
                        nr_self_collection[i, k] - nr_accretion_graupel[i, k] - ng_rain_snow[i, k]) * precip_fraction[i, k]
        qr_out[i, k] = qr[i, k]
        nr_per_volume[i, k] = nr[i, k] * ρ[i, k]
        qs_out[i, k] = qs[i, k]
        ns_per_volume[i, k] = ns[i, k] * ρ[i, k]
        qg_out[i, k] = qg[i, k]
        ng_per_volume[i, k] = ng[i, k] * ρ[i, k]
    end

    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qr_precip[i, k], nr_precip[i, k])
        λr[i, k], N₀r[i, k], nr_precip[i, k] = rain.slope, rain.intercept, rain.number
    end
    for k in 1:levels, i in 1:columns
        radius_rain_and_cloud[i, k] = add_rain_cloud_effective_radius(radius_rain_and_cloud[i, k], λr[i, k], N₀r[i, k],
            λc[i, k], μc[i, k], qr_precip[i, k], qc_cloud[i, k], nc_cloud[i, k])
    end
    return nothing
end
