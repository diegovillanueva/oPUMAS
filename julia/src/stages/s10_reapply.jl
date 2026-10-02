# ♻️🔁 Stage 10, re-apply the instantaneous adjustments (micro_pumas_v1.F90 lines
#    3172–3220). Stage 2 and 3 changed the state before the process rates; here the state
#    returns to its start-of-step values and the changes move into the tendencies.

"""
    reapply_instant_changes!(work, parameters)

♻️ Stage 10. Droplet activation, snow and graupel melting and rain freezing become
tendencies; snow sublimation and production join `nevapr` and `prain`.
"""
function reapply_instant_changes!(work, parameters)
    (; inputs, state, process_rates, tendencies, diagnostics, scratch) = work
    qs_sublimation_mean = process_rates.snow.sublimation.mass  # evapsnow
    qs_production = process_rates.snow.production.mass  # prodsnow
    nc = state.cloud_liquid.N
    qs = state.snow.q
    ns = state.snow.N
    qr = state.rain.q
    nr = state.rain.N
    qg = state.graupel.q
    ng = state.graupel.N
    nc_in = inputs.cloud_liquid.N  # ncn
    qs_in = inputs.snow.q  # qsn
    ns_in = inputs.snow.N  # nsn
    qr_in = inputs.rain.q  # qrn
    nr_in = inputs.rain.N  # nrn
    qg_in = inputs.graupel.q  # qgr
    ng_in = inputs.graupel.N  # ngr
    droplet_activation = inputs.activated_droplets  # npccn
    before = scratch.number_before  # dum_2D, here the value before the reset
    nc_tendency = tendencies.cloud_liquid.N  # nctend
    qs_tendency = tendencies.snow.q  # qstend
    ns_tendency = tendencies.snow.N  # nstend
    qr_tendency = tendencies.rain.q  # qrtend
    nr_tendency = tendencies.rain.N  # nrtend
    qg_tendency = tendencies.graupel.q  # qgtend
    ng_tendency = tendencies.graupel.N  # ngtend
    precip_evaporation = diagnostics.precipitation.evaporation  # nevapr
    precip_production = diagnostics.precipitation.production  # prain
    inv_Δt = 1.0 / work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        nc[i, k] = nc_in[i, k]
        nc_tendency[i, k] = nc_tendency[i, k] + droplet_activation[i, k]
        before[i, k] = qs[i, k]
        qs[i, k] = qs_in[i, k]
        qs_tendency[i, k] = qs_tendency[i, k] + (before[i, k] - qs[i, k]) * inv_Δt
        before[i, k] = ns[i, k]
        ns[i, k] = ns_in[i, k]
        ns_tendency[i, k] = ns_tendency[i, k] + (before[i, k] - ns[i, k]) * inv_Δt
        before[i, k] = qr[i, k]
        qr[i, k] = qr_in[i, k]
        qr_tendency[i, k] = qr_tendency[i, k] + (before[i, k] - qr[i, k]) * inv_Δt
        before[i, k] = nr[i, k]
        nr[i, k] = nr_in[i, k]
        nr_tendency[i, k] = nr_tendency[i, k] + (before[i, k] - nr[i, k]) * inv_Δt
        before[i, k] = qg[i, k]
        qg[i, k] = qg_in[i, k]
        qg_tendency[i, k] = qg_tendency[i, k] + (before[i, k] - qg[i, k]) * inv_Δt
        before[i, k] = ng[i, k]
        ng[i, k] = ng_in[i, k]
        ng_tendency[i, k] = ng_tendency[i, k] + (before[i, k] - ng[i, k]) * inv_Δt
        precip_evaporation[i, k] = precip_evaporation[i, k] + qs_sublimation_mean[i, k]
        precip_production[i, k] = precip_production[i, k] + qs_production[i, k]
    end
    return nothing
end
