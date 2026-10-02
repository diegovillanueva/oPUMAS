# ⚖️🛡️ Stage 8, conservation limiters (micro_pumas_v1.F90 lines 2584–2972). Every sink of
#    a species is scaled by one ratio so the species cannot go negative in the time step;
#    then evaporation and sublimation are limited so they do not supersaturate the air.

"""
    conservation_limiters!(work, parameters)

⚖️ Stage 8. One function per limiter, called in the Fortran order. Every limiter acts
point by point, so running them as separate loops gives the same numbers as the Fortran's
fused loops.
"""
function conservation_limiters!(work, parameters)
    (; settings) = parameters
    r = limiter_rates(work.rates)
    graupel_on = settings.do_hail || settings.do_graupel
    if settings.ifs_evaporation_scaling
        @. r.qr_evaporation = 0.15 * r.qr_evaporation
    end
    limit_cloud_liquid_mass!(work, r)
    limit_vapor_deposition!(work, parameters, r)
    limit_cloud_liquid_number!(work, r)
    small_rain_freezes_to_cloud_ice!(work, parameters, r)
    limit_rain_mass!(work, r)
    limit_rain_number!(work, r)
    if settings.prognostic_cloud_ice
        limit_cloud_ice_mass!(work, r)
        limit_cloud_ice_number!(work, settings, r)
    end
    limit_snow_mass!(work, r, graupel_on)
    limit_snow_number!(work, r, graupel_on)
    graupel_on && limit_graupel_mass!(work, r)
    limit_evaporation_to_saturation!(work, parameters, r)
    return nothing
end

"🗂️ The local rates the limiters read and rescale, under the names of the stages."
limiter_rates(rates) = (
    qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass,  # prc
    nc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.number,  # nprc1
    nr_autoconversion = rates.rain.autoconversion_from_cloud_liquid.number,  # nprc
    qc_accretion = rates.cloud_liquid.accretion_by_rain.mass,  # pra
    nc_accretion = rates.cloud_liquid.accretion_by_rain.number,  # npra
    qc_riming_snow = rates.cloud_liquid.riming_by_snow.mass,  # psacws
    nc_riming_snow = rates.cloud_liquid.riming_by_snow.number,  # npsacws
    qc_riming_graupel = rates.cloud_liquid.riming_by_graupel.mass,  # psacwg
    nc_riming_graupel = rates.cloud_liquid.riming_by_graupel.number,  # npsacwg
    qc_immersion_freezing = rates.cloud_liquid.immersion_freezing.mass,  # mnuccc
    nc_immersion_freezing = rates.cloud_liquid.immersion_freezing.number,  # nnuccc
    qc_contact_freezing = rates.cloud_liquid.contact_freezing.mass,  # mnucct
    nc_contact_freezing = rates.cloud_liquid.contact_freezing.number,  # nnucct
    qc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.mass,  # mnudep
    nc_deposition_nucleation = rates.cloud_liquid.deposition_nucleation.number,  # nnudep
    qc_bergeron_ice = rates.cloud_liquid.bergeron_to_cloud_ice.mass,  # berg
    qc_bergeron_snow = rates.cloud_liquid.bergeron_to_snow.mass,  # bergs
    nc_evaporation = rates.cloud_liquid.evaporation.number,  # nsubc
    nr_self_collection = rates.rain.self_collection.number,  # nragg
    qr_evaporation = rates.rain.evaporation.mass,  # pre
    nr_evaporation = rates.rain.evaporation.number,  # nsubr
    qr_accretion_snow = rates.rain.accretion_by_snow.mass,  # pracs
    nr_accretion_snow = rates.rain.accretion_by_snow.number,  # npracs
    qr_accretion_graupel = rates.rain.accretion_by_graupel.mass,  # pracg
    nr_accretion_graupel = rates.rain.accretion_by_graupel.number,  # npracg
    qr_freezing_snow = rates.rain.freezing_to_snow.mass,  # mnuccr
    nr_freezing_snow = rates.rain.freezing_to_snow.number,  # nnuccr
    qr_freezing_ice = rates.rain.freezing_to_cloud_ice.mass,  # mnuccri
    nr_freezing_ice = rates.rain.freezing_to_cloud_ice.number,  # nnuccri
    qi_nucleation = rates.cloud_ice.nucleation_from_vapor.mass,  # mnuccd
    ni_nucleation = rates.cloud_ice.nucleation_from_vapor.number,  # nnuccd
    qi_deposition = rates.cloud_ice.deposition.mass,  # vap_dep
    qi_sublimation = rates.cloud_ice.sublimation.mass,  # ice_sublim
    ni_sublimation = rates.cloud_ice.sublimation.number,  # nsubi
    qi_autoconversion = rates.cloud_ice.autoconversion_to_snow.mass,  # prci
    ni_autoconversion = rates.cloud_ice.autoconversion_to_snow.number,  # nprci
    qi_accretion_snow = rates.cloud_ice.accretion_by_snow.mass,  # prai
    ni_accretion_snow = rates.cloud_ice.accretion_by_snow.number,  # nprai
    qi_splintering_snow = rates.cloud_ice.splintering_by_snow.mass,  # msacwi
    ni_splintering_snow = rates.cloud_ice.splintering_by_snow.number,  # nsacwi
    qi_splintering_graupel = rates.cloud_ice.splintering_by_graupel.mass,  # qmultg
    ni_splintering_graupel = rates.cloud_ice.splintering_by_graupel.number,  # nmultg
    qi_splintering_rain = rates.cloud_ice.splintering_rain_graupel.mass,  # qmultrg
    ni_splintering_rain = rates.cloud_ice.splintering_rain_graupel.number,  # nmultrg
    qs_deposition = rates.snow.deposition.mass,  # vap_deps
    qs_sublimation = rates.snow.sublimation.mass,  # prds
    ns_sublimation = rates.snow.sublimation.number,  # nsubs
    ns_aggregation = rates.snow.aggregation.number,  # nsagg
    qs_collection_rain = rates.snow.collection_by_rain.mass,  # psacr
    qg_rimed_snow = rates.graupel.from_snow_riming.mass,  # pgsacw
    ng_rimed_snow = rates.graupel.from_snow_riming.number,  # nscng
    qg_rain_snow = rates.graupel.rain_snow_collection.mass,  # pgracs
    ng_rain_snow = rates.graupel.rain_snow_collection.number,  # ngracs
    qg_deposition = rates.graupel.deposition.mass  # prdg
)

"✖️ Multiply every rate at point `(i, k)` by the limiter `ratio`."
function scale!(ratio, i, k, rates...)
    for rate in rates
        rate[i, k] = rate[i, k] * ratio
    end
    return nothing
end

"⏱️🔁 Time step and its inverse (`deltat`, `rdeltat`)."
time_step(work) = (work.inputs.time_step, 1.0 / work.inputs.time_step)

"""
💧 Cloud liquid mass: if its sinks would remove more than `qc` in one step, all of them
are scaled by one ratio, kept as the diagnostic `qcrat`.
"""
function limit_cloud_liquid_mass!(work, r)
    (; qc_autoconversion, qc_accretion, qc_immersion_freezing, qc_contact_freezing,
    qi_splintering_snow, qc_riming_snow, qc_bergeron_snow, qi_splintering_graupel,
    qc_riming_graupel, qg_rimed_snow, qc_bergeron_ice) = r
    qc = work.state.cloud_liquid.q
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    qc_limiter_ratio = work.diagnostics.cloud_liquid.limiter_ratio  # qcrat
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        sink_rate = (qc_autoconversion[i, k] + qc_accretion[i, k] +
                     qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                     qi_splintering_snow[i, k] +
                     qc_riming_snow[i, k] + qc_bergeron_snow[i, k] +
                     qi_splintering_graupel[i, k] + qc_riming_graupel[i, k] +
                     qg_rimed_snow[i, k]) * liquid_fraction[i, k] + qc_bergeron_ice[i, k]
        if sink_rate * Δt > qc[i, k]
            ratio = qc[i, k] * inv_Δt / sink_rate * one_minus_small
            scale!(ratio, i, k, qi_splintering_graupel, qc_riming_graupel, qg_rimed_snow,
                qc_autoconversion, qc_accretion, qc_immersion_freezing,
                qc_contact_freezing, qi_splintering_snow, qc_riming_snow,
                qc_bergeron_snow, qc_bergeron_ice)
            qc_limiter_ratio[i, k] = ratio
        else
            qc_limiter_ratio[i, k] = 1.0
        end
    end
    return nothing
end

"""
🧊🌫️ Vapor deposition onto ice and snow: reduced where the cloud liquid limiter acted,
then capped so it cannot take the air below ice saturation in one step.
"""
function limit_vapor_deposition!(work, parameters, r)
    (; qi_deposition, qi_nucleation, qs_deposition) = r
    (; cp, Rv, Ls) = host_constants(parameters)
    Ls_squared = Ls^2
    T = work.inputs.T  # t
    qv = work.inputs.q_vapor  # q
    qsat_ice = work.state.air.q_sat_ice  # qvi
    qc = work.state.cloud_liquid.q
    qc_limiter_ratio = work.diagnostics.cloud_liquid.limiter_ratio  # qcrat
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= smallest_mixing_ratio
            qi_deposition[i, k] = qi_deposition[i, k] * (1.0 - qc_limiter_ratio[i, k])
            qs_deposition[i, k] = qs_deposition[i, k] * (1.0 - qc_limiter_ratio[i, k])
        end
        qi_nucleation[i, k] = max(0.0, qi_nucleation[i, k])
        qi_deposition[i, k] = max(0.0, qi_deposition[i, k])
        qs_deposition[i, k] = max(0.0, qs_deposition[i, k])
        deposition_total = qi_deposition[i, k] + qi_nucleation[i, k] + qs_deposition[i, k]
        deposition_total > 1.0e-20 || continue
        deposition_max = (qv[i, k] - qsat_ice[i, k]) /
                         (1.0 +
                          Ls_squared * qsat_ice[i, k] / (cp * Rv * (T[i, k] * T[i, k]))) *
                         inv_Δt
        deposition_max = max(deposition_max, 0.0)
        if deposition_total > deposition_max
            qi_nucleation[i, k] = deposition_max * qi_nucleation[i, k] / deposition_total
            qi_deposition[i, k] = deposition_max * qi_deposition[i, k] / deposition_total
            qs_deposition[i, k] = deposition_max * qs_deposition[i, k] / deposition_total
        end
    end
    return nothing
end

"💧🔢 Cloud droplet number: all sinks scaled by one ratio if they exceed `nc`."
function limit_cloud_liquid_number!(work, r)
    (; nc_autoconversion, nc_accretion, nc_immersion_freezing, nc_contact_freezing,
    nc_riming_snow, nc_evaporation, nc_riming_graupel) = r
    nc = work.state.cloud_liquid.N
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        sink_rate = (nc_autoconversion[i, k] + nc_accretion[i, k] +
                     nc_immersion_freezing[i, k] + nc_contact_freezing[i, k] +
                     nc_riming_snow[i, k] - nc_evaporation[i, k] +
                     nc_riming_graupel[i, k]) * liquid_fraction[i, k]
        if sink_rate * Δt > nc[i, k]
            ratio = nc[i, k] * inv_Δt / sink_rate * one_minus_small
            scale!(ratio, i, k, nc_riming_graupel, nc_autoconversion, nc_accretion,
                nc_immersion_freezing, nc_contact_freezing, nc_riming_snow,
                nc_evaporation)
        end
    end
    return nothing
end

"""
🌧️🧊 Rain freezing goes to cloud ice instead of snow when the drops are smaller than the
ice to snow autoconversion size (`1/λr < dcs`, prognostic ice only).
"""
function small_rain_freezes_to_cloud_ice!(work, parameters, r)
    (; qr_freezing_ice, nr_freezing_ice, qr_freezing_snow, nr_freezing_snow) = r
    (; settings) = parameters
    λr = work.incloud.rain.λ  # lamr
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
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
    return nothing
end

"""
🌧️ Rain mass: the sinks are scaled if they exceed `qr` plus what autoconversion and
accretion add in the step. Sink sums below `qsmall` are never limited (FINDINGS.md).
"""
function limit_rain_mass!(work, r)
    (; qr_evaporation, qr_accretion_snow, qr_freezing_snow, qr_freezing_ice,
    qi_splintering_rain, qr_accretion_graupel, qg_rain_snow, qc_accretion,
    qc_autoconversion) = r
    qr = work.state.rain.q
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        sinks = -qr_evaporation[i, k] + qr_accretion_snow[i, k] + qr_freezing_snow[i, k] +
                qr_freezing_ice[i, k] + qi_splintering_rain[i, k] +
                qr_accretion_graupel[i, k] + qg_rain_snow[i, k]
        from_cloud = (qc_accretion[i, k] + qc_autoconversion[i, k]) * liquid_fraction[i, k]
        sink_over_step = (sinks * precip_fraction[i, k] - from_cloud) * Δt
        if sink_over_step > qr[i, k] && sinks >= smallest_mixing_ratio
            ratio = (qr[i, k] * inv_Δt + from_cloud) / precip_fraction[i, k] /
                    sinks * one_minus_small
            scale!(ratio, i, k, qi_splintering_rain, qr_accretion_graupel, qg_rain_snow,
                qr_evaporation, qr_accretion_snow, qr_freezing_snow, qr_freezing_ice)
        end
    end
    return nothing
end

"""
🌧️🔢 Rain number: evaporation removes drops in proportion to mass, then the sinks are
scaled if they exceed `nr` plus the drops autoconversion adds.
"""
function limit_rain_number!(work, r)
    (; qr_evaporation, nr_evaporation, nr_accretion_snow, nr_freezing_snow,
    nr_freezing_ice, nr_self_collection, nr_accretion_graupel, ng_rain_snow,
    nr_autoconversion) = r
    qr = work.state.rain.q
    nr = work.state.rain.N
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
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
            ratio = nr_sources / precip_fraction[i, k] / nr_sinks * one_minus_small
            scale!(ratio, i, k, nr_accretion_graupel, ng_rain_snow, nr_self_collection,
                nr_accretion_snow, nr_freezing_snow, nr_evaporation, nr_freezing_ice)
        end
    end
    return nothing
end

"""
🧊 Cloud ice mass: only autoconversion to snow, collection by snow and sublimation are
sinks; they are scaled if they exceed `qi` plus all sources of the step.
"""
function limit_cloud_ice_mass!(work, r)
    (; qc_immersion_freezing, qc_contact_freezing, qc_deposition_nucleation,
    qi_splintering_snow, qi_splintering_graupel, qi_autoconversion, qi_accretion_snow,
    qi_splintering_rain, qr_freezing_ice, qi_sublimation, qi_deposition,
    qc_bergeron_ice, qi_nucleation) = r
    qi = work.state.cloud_ice.q
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    ice_fraction = work.incloud.ice_fraction  # icldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        from_liquid = (qc_immersion_freezing[i, k] + qc_contact_freezing[i, k] +
                       qc_deposition_nucleation[i, k] + qi_splintering_snow[i, k] +
                       qi_splintering_graupel[i, k]) * liquid_fraction[i, k]
        from_rain = (qi_splintering_rain[i, k] + qr_freezing_ice[i, k]) *
                    precip_fraction[i, k]
        to_snow = (qi_autoconversion[i, k] + qi_accretion_snow[i, k]) * ice_fraction[i, k]
        sink_over_step = (-from_liquid + to_snow - from_rain - qi_sublimation[i, k] -
                          qi_deposition[i, k] - qc_bergeron_ice[i, k] -
                          qi_nucleation[i, k]) * Δt
        if sink_over_step > qi[i, k]
            ratio = (qi[i, k] * inv_Δt + qi_deposition[i, k] + qc_bergeron_ice[i, k] +
                     qi_nucleation[i, k] + from_liquid + from_rain) /
                    (to_snow - qi_sublimation[i, k]) * one_minus_small
            scale!(ratio, i, k, qi_autoconversion, qi_accretion_snow, qi_sublimation)
        end
    end
    return nothing
end

"""
🧊🔢 Cloud ice number, like the mass. Host immersion freezing counts as a source only
with `host_heterogeneous_freezing` (`use_hetfrz_classnuc`).
"""
function limit_cloud_ice_number!(work, settings, r)
    (; nc_contact_freezing, nc_immersion_freezing, nc_deposition_nucleation,
    ni_splintering_snow, ni_splintering_graupel, ni_autoconversion, ni_accretion_snow,
    ni_sublimation, ni_splintering_rain, nr_freezing_ice, ni_nucleation) = r
    ni = work.state.cloud_ice.N
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    ice_fraction = work.incloud.ice_fraction  # icldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        nc_immersion_counted = settings.host_heterogeneous_freezing ?
                               nc_immersion_freezing[i, k] : 0.0
        from_liquid = (nc_contact_freezing[i, k] + nc_immersion_counted +
                       nc_deposition_nucleation[i, k] + ni_splintering_snow[i, k] +
                       ni_splintering_graupel[i, k]) * liquid_fraction[i, k]
        from_rain = (nr_freezing_ice[i, k] + ni_splintering_rain[i, k]) *
                    precip_fraction[i, k]
        sinks = (ni_autoconversion[i, k] + ni_accretion_snow[i, k] -
                 ni_sublimation[i, k]) * ice_fraction[i, k]
        sink_over_step = (-from_liquid + sinks - from_rain - ni_nucleation[i, k]) * Δt
        if sink_over_step > ni[i, k]
            ratio = (ni[i, k] * inv_Δt + ni_nucleation[i, k] + from_liquid + from_rain) /
                    sinks * one_minus_small
            scale!(ratio, i, k, ni_autoconversion, ni_accretion_snow, ni_sublimation)
        end
    end
    return nothing
end

"""
❄️ Snow mass: only sublimation (and, with graupel, collection by rain) is limited. Without
graupel, rain freezing feeds snow; with graupel it feeds graupel.
"""
function limit_snow_mass!(work, r, graupel_on)
    (; qs_sublimation, qr_accretion_snow, qs_collection_rain, qr_freezing_snow,
    qi_accretion_snow, qi_autoconversion, qc_bergeron_snow, qc_riming_snow,
    qs_deposition) = r
    qs = work.state.snow.q
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    ice_fraction = work.incloud.ice_fraction  # icldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        pf = precip_fraction[i, k]
        from_ice = (qi_accretion_snow[i, k] + qi_autoconversion[i, k]) * ice_fraction[i, k]
        from_liquid = (qc_bergeron_snow[i, k] + qc_riming_snow[i, k]) *
                      liquid_fraction[i, k]
        if graupel_on
            sink_over_step = (-(qs_sublimation[i, k] + qr_accretion_snow[i, k] -
                                qs_collection_rain[i, k]) * pf -
                              from_ice - from_liquid - qs_deposition[i, k]) * Δt
        else
            sink_over_step = (-(qs_sublimation[i, k] + qr_accretion_snow[i, k] +
                                qr_freezing_snow[i, k]) * pf -
                              from_ice - from_liquid - qs_deposition[i, k]) * Δt
        end
        sink_over_step > qs[i, k] || continue
        (qs_collection_rain[i, k] - qs_sublimation[i, k]) >= smallest_mixing_ratio ||
            continue
        if graupel_on
            ratio = (qs[i, k] * inv_Δt + from_ice + from_liquid + qs_deposition[i, k] +
                     qr_accretion_snow[i, k] * pf) / pf /
                    (qs_collection_rain[i, k] - qs_sublimation[i, k]) * one_minus_small
            qs_collection_rain[i, k] = qs_collection_rain[i, k] * ratio
        else
            ratio = (qs[i, k] * inv_Δt + from_ice + from_liquid + qs_deposition[i, k] +
                     (qr_accretion_snow[i, k] + qr_freezing_snow[i, k]) * pf) / pf /
                    (-qs_sublimation[i, k]) * one_minus_small
        end
        qs_sublimation[i, k] = qs_sublimation[i, k] * ratio
    end
    return nothing
end

"❄️🔢 Snow number: sublimation removes no flakes; aggregation and graupel formation do."
function limit_snow_number!(work, r, graupel_on)
    (; ns_sublimation, ns_aggregation, ng_rain_snow, ng_rimed_snow, nr_freezing_snow,
    ni_autoconversion) = r
    ns = work.state.snow.N
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    ice_fraction = work.incloud.ice_fraction  # icldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        pf = precip_fraction[i, k]
        ns_sublimation[i, k] = 0.0
        from_ice = ni_autoconversion[i, k] * ice_fraction[i, k]
        if graupel_on
            sink_over_step = ((-ns_aggregation[i, k] - ns_sublimation[i, k] +
                               ng_rain_snow[i, k]) * pf -
                              from_ice + ng_rimed_snow[i, k] * liquid_fraction[i, k]) * Δt
        else
            sink_over_step = ((-ns_aggregation[i, k] - ns_sublimation[i, k] -
                               nr_freezing_snow[i, k]) * pf - from_ice) * Δt
        end
        sink_over_step > ns[i, k] || continue
        if graupel_on
            ratio = (ns[i, k] * inv_Δt + from_ice) / pf /
                    (-ns_sublimation[i, k] - ns_aggregation[i, k] + ng_rain_snow[i, k] +
                     liquid_fraction[i, k] / pf * ng_rimed_snow[i, k]) * one_minus_small
            scale!(ratio, i, k, ng_rimed_snow, ng_rain_snow)
        else
            ratio = (ns[i, k] * inv_Δt + nr_freezing_snow[i, k] * pf + from_ice) / pf /
                    (-ns_sublimation[i, k] - ns_aggregation[i, k]) * one_minus_small
        end
        scale!(ratio, i, k, ns_sublimation, ns_aggregation)
    end
    return nothing
end

"🟤 Graupel mass: only its sublimation is limited (graupel number has no sinks)."
function limit_graupel_mass!(work, r)
    (; qr_accretion_graupel, qg_rain_snow, qg_deposition, qs_collection_rain,
    qr_freezing_snow, qc_riming_graupel, qg_rimed_snow) = r
    qg = work.state.graupel.q
    liquid_fraction = work.incloud.liquid_fraction  # lcldm
    precip_fraction = work.incloud.precipitation_fraction  # precip_frac
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        pf = precip_fraction[i, k]
        sink_over_step = ((-qr_accretion_graupel[i, k] - qg_rain_snow[i, k] -
                           qg_deposition[i, k] - qs_collection_rain[i, k] -
                           qr_freezing_snow[i, k]) * pf +
                          (-qc_riming_graupel[i, k] - qg_rimed_snow[i, k]) *
                          liquid_fraction[i, k]) * Δt
        sink_over_step > qg[i, k] || continue
        ratio = (qg[i, k] * inv_Δt +
                 (qr_accretion_graupel[i, k] + qg_rain_snow[i, k] +
                  qs_collection_rain[i, k] + qr_freezing_snow[i, k]) * pf +
                 (qc_riming_graupel[i, k] + qg_rimed_snow[i, k]) * liquid_fraction[i, k]) /
                ((-qg_deposition[i, k]) * pf) * one_minus_small
        qg_deposition[i, k] = qg_deposition[i, k] * ratio
    end
    return nothing
end

"""
🌫️ Evaporation and sublimation may not take the air past saturation. Where they would,
their total is reset to the vapor deficit after deposition, over liquid for rain and over
ice for snow, graupel and cloud ice, shared out in the old proportions.
"""
function limit_evaporation_to_saturation!(work, parameters, r)
    (; qr_evaporation, qs_sublimation, qg_deposition, qi_sublimation, qi_deposition,
    qi_nucleation, qs_deposition) = r
    (; saturation) = parameters
    (; cp, Rv, Lv, Ls) = host_constants(parameters)
    Lv_squared, Ls_squared = Lv^2, Ls^2
    (; inputs, incloud, scratch) = work
    T = inputs.T  # t
    qv = inputs.q_vapor  # q
    p = inputs.pressure
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    qv_after_deposition = scratch.number_before  # dum_2D, here the vapor after deposition
    T_estimate = scratch.temperature_estimate  # ttmpA
    qv_estimate = scratch.vapor_estimate  # qtmpAI
    qsat_estimate = scratch.q_sat_estimate  # qvnAI
    qsat_adjusted = scratch.q_sat_adjusted  # qvnA
    rain_share = scratch.evaporation_share.rain  # dum1A
    snow_share = scratch.evaporation_share.snow  # dum2A
    graupel_share = scratch.evaporation_share.graupel  # dum3A
    Δt, inv_Δt = time_step(work)
    columns, levels = grid_size(work)

    # 🔮 vapor and temperature if every rate ran unlimited
    for k in 1:levels, i in 1:columns
        qv_estimate[i, k] = qv[i, k] -
                            (qi_sublimation[i, k] + qi_deposition[i, k] +
                             qi_nucleation[i, k] + qs_deposition[i, k] +
                             (qr_evaporation[i, k] + qs_sublimation[i, k] +
                              qg_deposition[i, k]) * precip_fraction[i, k]) * Δt
        T_estimate[i, k] = T[i, k] +
                           ((qr_evaporation[i, k] * precip_fraction[i, k]) * Lv +
                            ((qs_sublimation[i, k] + qg_deposition[i, k]) *
                             precip_fraction[i, k] +
                             qi_deposition[i, k] + qs_deposition[i, k] +
                             qi_sublimation[i, k] + qi_nucleation[i, k]) * Ls) * Δt / cp
    end
    for k in 1:levels, i in 1:columns
        qsat_estimate[i, k] = saturation_over_water(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    function evaporation(i, k)
        (qr_evaporation[i, k] + qs_sublimation[i, k] +
         qg_deposition[i, k]) * precip_fraction[i, k] + qi_sublimation[i, k]
    end
    function supersaturating(i, k)
        evaporation(i, k) < -1.0e-20 &&
            qv_estimate[i, k] > qsat_estimate[i, k]
    end

    # ⚖️ shares of the evaporation, and the state after deposition alone
    for k in 1:levels, i in 1:columns
        supersaturating(i, k) || continue
        total = evaporation(i, k)
        rain_share[i, k] = qr_evaporation[i, k] * precip_fraction[i, k] / total
        snow_share[i, k] = qs_sublimation[i, k] * precip_fraction[i, k] / total
        graupel_share[i, k] = qg_deposition[i, k] * precip_fraction[i, k] / total
        T_estimate[i, k] = T[i, k] +
                           ((qi_deposition[i, k] + qs_deposition[i, k] +
                             qi_nucleation[i, k]) * Ls) * Δt / cp
        qv_after_deposition[i, k] = qv[i, k] -
                                    (qi_deposition[i, k] + qs_deposition[i, k] +
                                     qi_nucleation[i, k]) * Δt
    end

    # 💧 rain evaporates at most the deficit to liquid saturation
    for k in 1:levels, i in 1:columns
        qsat_adjusted[i, k] = saturation_over_water(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        supersaturating(i, k) || continue
        vapor_excess = (qv_after_deposition[i, k] - qsat_adjusted[i, k]) /
                       (1.0 +
                        Lv_squared * qsat_adjusted[i, k] /
                        (cp * Rv * (T_estimate[i, k] * T_estimate[i, k])))
        vapor_excess = min(vapor_excess, 0.0)
        qr_evaporation[i, k] = vapor_excess * rain_share[i, k] * inv_Δt /
                               precip_fraction[i, k]
    end

    # 🧊 snow, graupel and cloud ice sublimate at most the deficit to ice saturation
    for k in 1:levels, i in 1:columns
        qsat_adjusted[i, k] = saturation_over_ice(
            saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end
    for k in 1:levels, i in 1:columns
        supersaturating(i, k) || continue
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
    return nothing
end
