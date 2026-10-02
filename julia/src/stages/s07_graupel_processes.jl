# 🟤🌧️ Stage 7, process rates of graupel or hail and evaporation of precipitation
#    (micro_pumas_v1.F90 lines 2531–2583).

"""
    graupel_processes!(work, parameters)

🟤 Stage 7. Collection of snow, cloud liquid and rain by graupel, rime splintering, and
evaporation and sublimation of rain, snow and graupel.
"""
function graupel_processes!(work, parameters)
    (; settings) = parameters
    if !(settings.do_hail || settings.do_graupel)
        error("precipitation evaporation without graupel or hail is not translated yet")
    end
    dense = dense_precipitation(parameters)
    (; inputs, state, incloud, rates, process_rates, diagnostics) = work
    T = inputs.T  # t
    ρ = state.air.density  # rho
    a_rain = state.air.fall_speed_prefactor.rain  # arn
    a_snow = state.air.fall_speed_prefactor.snow  # asn
    a_graupel = state.air.fall_speed_prefactor.graupel  # agn
    liquid_fraction = incloud.liquid_fraction  # lcldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    qc_cloud = incloud.cloud_liquid.q  # qcic
    nc_cloud = incloud.cloud_liquid.N  # ncic
    qi_cloud = incloud.cloud_ice.q  # qiic
    qr_precip = incloud.rain.q  # qric
    nr_precip = incloud.rain.N  # nric
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    vr_number = incloud.rain.fall_speed_number  # unr
    vr_mass = process_rates.rain.fall_speed.mass  # umr
    qs_precip = incloud.snow.q  # qsic
    ns_precip = incloud.snow.N  # nsic
    λs = incloud.snow.λ  # lams
    N₀s = incloud.snow.N₀  # n0s
    vs_mass = process_rates.snow.fall_speed.mass  # ums
    qg_precip = incloud.graupel.q  # qgic
    λg = incloud.graupel.λ  # lamg
    N₀g = incloud.graupel.N₀  # n0g
    vg_number = incloud.graupel.fall_speed_number  # ung
    vg_mass = process_rates.graupel.fall_speed.mass  # umg
    qs_collection_rain = rates.snow.collection_by_rain.mass  # psacr
    qc_riming_graupel = rates.cloud_liquid.riming_by_graupel.mass  # psacwg
    nc_riming_graupel = rates.cloud_liquid.riming_by_graupel.number  # npsacwg
    qc_riming_snow = rates.cloud_liquid.riming_by_snow.mass  # psacws
    qg_rimed_snow = rates.graupel.from_snow_riming.mass  # pgsacw
    ng_rimed_snow = rates.graupel.from_snow_riming.number  # nscng
    qr_accretion_graupel = rates.rain.accretion_by_graupel.mass  # pracg
    nr_accretion_graupel = rates.rain.accretion_by_graupel.number  # npracg
    qr_accretion_snow = rates.rain.accretion_by_snow.mass  # pracs
    nr_accretion_snow = rates.rain.accretion_by_snow.number  # npracs
    qg_rain_snow = rates.graupel.rain_snow_collection.mass  # pgracs
    ng_rain_snow = rates.graupel.rain_snow_collection.number  # ngracs
    qi_splintering_graupel = rates.cloud_ice.splintering_by_graupel.mass  # qmultg
    ni_splintering_graupel = rates.cloud_ice.splintering_by_graupel.number  # nmultg
    qi_splintering_rain = rates.cloud_ice.splintering_rain_graupel.mass  # qmultrg
    ni_splintering_rain = rates.cloud_ice.splintering_rain_graupel.number  # nmultrg
    qr_evaporation = rates.rain.evaporation.mass  # pre
    qs_sublimation = rates.snow.sublimation.mass  # prds
    qg_deposition = rates.graupel.deposition.mass  # prdg
    evaporation_area = diagnostics.precipitation.evaporation_area  # am_evp_st
    Δt = work.inputs.time_step
    iaccr = settings.ice_accretion_factor
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        qs_collection_rain[i, k] = snow_collection_by_rain(
            vr_mass[i, k], vs_mass[i, k], ρ[i, k], qs_precip[i, k],
            qr_precip[i, k], λr[i, k], N₀r[i, k], λs[i, k], N₀s[i, k])
    end
    for k in 1:levels, i in 1:columns
        rate = droplet_collection_by_graupel(
            qg_precip[i, k], qc_cloud[i, k], nc_cloud[i, k], ρ[i, k],
            N₀g[i, k], λg[i, k], dense.b, a_graupel[i, k])
        qc_riming_graupel[i, k], nc_riming_graupel[i, k] = rate.mass, rate.number
    end
    for k in 1:levels, i in 1:columns
        qc_riming_graupel[i, k] = qc_riming_graupel[i, k] * iaccr
        nc_riming_graupel[i, k] = nc_riming_graupel[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        rate = rimed_snow_to_graupel(
            parameters, qc_riming_snow[i, k], qs_precip[i, k], qc_cloud[i, k],
            ns_precip[i, k], ρ[i, k], ρ_snow, dense.density, a_snow[i, k], λs[i, k], N₀s[
                i, k],
            Δt)
        qg_rimed_snow[i, k], ng_rimed_snow[i, k], qc_riming_snow[i, k] = rate.mass,
        rate.number, rate.riming
    end
    for k in 1:levels, i in 1:columns
        rate = rain_collection_by_graupel(
            qr_precip[i, k], qg_precip[i, k], vg_mass[i, k], vr_mass[i, k],
            vg_number[i, k], vr_number[i, k], ρ[i, k], N₀r[i, k], λr[i, k], N₀g[i, k], λg[
                i, k])
        qr_accretion_graupel[i, k], nr_accretion_graupel[i, k] = rate.mass, rate.number
    end
    for k in 1:levels, i in 1:columns
        qr_accretion_graupel[i, k] = qr_accretion_graupel[i, k] * iaccr
        nr_accretion_graupel[i, k] = nr_accretion_graupel[i, k] * iaccr
    end
    for k in 1:levels, i in 1:columns
        rate = rain_snow_collisions_to_graupel(
            qr_accretion_snow[i, k], nr_accretion_snow[i, k], qs_collection_rain[i, k],
            qs_precip[i, k], qr_precip[i, k], nr_precip[i, k], ns_precip[i, k],
            N₀s[i, k], λs[i, k], N₀r[i, k], λr[i, k], Δt)
        qg_rain_snow[i, k], ng_rain_snow[i, k] = rate.mass, rate.number
        qr_accretion_snow[i, k], nr_accretion_snow[i, k] = rate.rain_by_snow,
        rate.rain_number_by_snow
        qs_collection_rain[i, k] = rate.snow_by_rain
    end
    for k in 1:levels, i in 1:columns
        rate = rime_splintering_by_graupel(
            parameters, T[i, k], qc_cloud[i, k], qr_precip[i, k],
            qg_precip[i, k], qc_riming_graupel[i, k], qr_accretion_graupel[i, k])
        qi_splintering_graupel[i, k], ni_splintering_graupel[i, k] = rate.droplet_mass,
        rate.droplet_number
        qi_splintering_rain[i, k], ni_splintering_rain[i, k] = rate.rain_mass,
        rate.rain_number
        qc_riming_graupel[i, k], qr_accretion_graupel[i, k] = rate.droplet_riming,
        rate.rain_riming
    end
    for k in 1:levels, i in 1:columns
        rain = (q = qr_precip[i, k], λ = λr[i, k], N₀ = N₀r[i, k], a = a_rain[i, k])
        snow = (q = qs_precip[i, k], λ = λs[i, k], N₀ = N₀s[i, k], a = a_snow[i, k])
        graupel = (q = qg_precip[i, k], λ = λg[i, k], N₀ = N₀g[i, k], a = a_graupel[i, k],
            b = dense.b)
        rate = precipitation_evaporation(
            parameters, air_at(work, i, k), liquid_fraction[i, k],
            precip_fraction[i, k], qc_cloud[i, k], qi_cloud[i, k], rain, snow, graupel)
        qr_evaporation[i, k], qs_sublimation[i, k], qg_deposition[i, k] = rate.rain,
        rate.snow, rate.graupel
        evaporation_area[i, k] = rate.area
    end
    return nothing
end
