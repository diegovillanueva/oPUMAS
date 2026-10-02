# 🌡️🔁 Stage 3, instantaneous melting and freezing, first pass (micro_pumas_v1.F90
#    lines 1784–1898). The state is changed directly; stage 10 moves it into tendencies.

"""
    instant_phase_change!(work, parameters)

🌡️ Stage 3. Snow and graupel melt above +2 °C and rain freezes below the rain freezing
temperature, each limited so that the latent heat does not move the temperature across
the threshold.
"""
function instant_phase_change!(work, parameters)
    (; settings) = parameters
    (; cp, Lf, T_snow_melt, T_rain_freeze) = host_constants(parameters)
    (; inputs, state, rates, process_rates, tendencies) = work
    T = inputs.T  # t
    qr = state.rain.q
    nr = state.rain.N
    qs = state.snow.q
    ns = state.snow.N
    qg = state.graupel.q
    ng = state.graupel.N
    heating = tendencies.heating  # tlat
    qs_instant_melting = rates.snow.instant_melting.mass  # minstsm
    ns_instant_melting = rates.snow.instant_melting.number  # ninstsm
    qg_instant_melting = rates.graupel.instant_melting.mass  # minstgm
    ng_instant_melting = rates.graupel.instant_melting.number  # ninstgm
    qr_instant_freezing = rates.rain.instant_freezing.mass  # minstrf
    nr_instant_freezing = rates.rain.instant_freezing.number  # ninstrf
    heating_snow_melting = process_rates.heating.snow_melting  # meltsdttot
    qs_melting_mean = process_rates.snow.melting.mass  # meltstot
    qg_melting_mean = process_rates.graupel.melting.mass  # meltgtot
    heating_rain_freezing = process_rates.heating.rain_freezing  # frzrdttot
    inv_Δt = 1.0 / work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        # ❄️➡️💧 melting of snow at +2 °C
        if T[i, k] > T_snow_melt && qs[i, k] > 0.0
            dum = -Lf / cp * qs[i, k]
            if T[i, k] + dum < T_snow_melt
                dum = (T[i, k] - T_snow_melt) * cp / Lf
                dum = dum / qs[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            qs_instant_melting[i, k] = dum * qs[i, k]
            ns_instant_melting[i, k] = dum * ns[i, k]
            dum1 = -Lf * qs_instant_melting[i, k] * inv_Δt
            heating[i, k] = heating[i, k] + dum1
            heating_snow_melting[i, k] = heating_snow_melting[i, k] + dum1
            qs_melting_mean[i, k] = qs_instant_melting[i, k] * inv_Δt
            qs[i, k] = max(qs[i, k] - qs_instant_melting[i, k], 0.0)
            ns[i, k] = max(ns[i, k] - ns_instant_melting[i, k], 0.0)
            qr[i, k] = max(qr[i, k] + qs_instant_melting[i, k], 0.0)
            nr[i, k] = max(nr[i, k] + ns_instant_melting[i, k], 0.0)
        end

        # 🟤➡️💧 melting of graupel at +2 °C
        if T[i, k] > T_snow_melt && qg[i, k] > 0.0
            dum = -Lf / cp * qg[i, k]
            if T[i, k] + dum < T_snow_melt
                dum = (T[i, k] - T_snow_melt) * cp / Lf
                dum = dum / qg[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            qg_instant_melting[i, k] = dum * qg[i, k]
            ng_instant_melting[i, k] = dum * ng[i, k]
            dum1 = -Lf * qg_instant_melting[i, k] * inv_Δt
            heating[i, k] = heating[i, k] + dum1
            heating_snow_melting[i, k] = heating_snow_melting[i, k] + dum1
            qg_melting_mean[i, k] = qg_instant_melting[i, k] * inv_Δt
            qg[i, k] = max(qg[i, k] - qg_instant_melting[i, k], 0.0)
            ng[i, k] = max(ng[i, k] - ng_instant_melting[i, k], 0.0)
            qr[i, k] = max(qr[i, k] + qg_instant_melting[i, k], 0.0)
            nr[i, k] = max(nr[i, k] + ng_instant_melting[i, k], 0.0)
        end

        # 💧➡️🧊 freezing of rain, to graupel (or hail) if on, else to snow
        if T[i, k] < T_rain_freeze && qr[i, k] > 0.0
            dum = Lf / cp * qr[i, k]
            if T[i, k] + dum > T_rain_freeze
                dum = -(T[i, k] - T_rain_freeze) * cp / Lf
                dum = dum / qr[i, k]
                dum = max(0.0, dum)
                dum = min(1.0, dum)
            else
                dum = 1.0
            end
            qr_instant_freezing[i, k] = dum * qr[i, k]
            nr_instant_freezing[i, k] = dum * nr[i, k]
            dum1 = Lf * qr_instant_freezing[i, k] * inv_Δt
            heating[i, k] = heating[i, k] + dum1
            heating_rain_freezing[i, k] = heating_rain_freezing[i, k] + dum1
            qr[i, k] = max(qr[i, k] - qr_instant_freezing[i, k], 0.0)
            nr[i, k] = max(nr[i, k] - nr_instant_freezing[i, k], 0.0)
            if settings.do_hail || settings.do_graupel
                qg[i, k] = max(qg[i, k] + qr_instant_freezing[i, k], 0.0)
                ng[i, k] = max(ng[i, k] + nr_instant_freezing[i, k], 0.0)
            else
                qs[i, k] = max(qs[i, k] + qr_instant_freezing[i, k], 0.0)
                ns[i, k] = max(ns[i, k] + nr_instant_freezing[i, k], 0.0)
            end
        end
    end
    return nothing
end
