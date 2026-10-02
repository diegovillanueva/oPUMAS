# 🌡️🧊 Stage 13, instantaneous processes after sedimentation (micro_pumas_v1.F90 lines
#    3733–4020): snow and graupel melting, rain freezing, cloud ice melting, homogeneous
#    droplet freezing, the ice number limit and the removal of supersaturation.

"""
    instant_phase_change_after_sedimentation!(work, parameters)

🌡️ Stage 13. Each instantaneous change is limited so that it does not push the
end-of-step temperature across its threshold; tendencies and `tlat` are updated.
"""
function instant_phase_change_after_sedimentation!(work, parameters)
    (; settings, saturation, hydrometeors) = parameters
    (; cp, Rv, Lv, Ls, Lf, T₀, T_snow_melt, T_rain_freeze) = host_constants(parameters)
    Lv_squared = Lv^2
    (; inputs, state, incloud, process_rates, tendencies, diagnostics, scratch) = work
    liquid_fraction = incloud.liquid_fraction  # lcldm
    heating_snow_melting = process_rates.heating.snow_melting  # meltsdttot
    heating_rain_freezing = process_rates.heating.rain_freezing  # frzrdttot
    qi_melting_mean = process_rates.cloud_ice.melting.mass  # melttot
    ni_melting_mean = process_rates.cloud_ice.melting.number  # nmelttot
    qc_homogeneous_freezing_mean = process_rates.cloud_liquid.homogeneous_freezing.mass  # homotot
    nc_homogeneous_freezing_mean = process_rates.cloud_liquid.homogeneous_freezing.number  # nhomotot
    qc_residual_condensation = process_rates.cloud_liquid.residual_condensation.mass  # qcrestot
    qi_residual_deposition = process_rates.cloud_ice.residual_deposition.mass  # qirestot
    qv_residual_condensation = process_rates.vapor.residual_condensation.mass  # qvres
    T = inputs.T  # t
    qv = inputs.q_vapor  # q
    p = inputs.pressure
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
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    heating = tendencies.heating  # tlat
    qv_tendency = tendencies.vapor  # qvlat
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
    net_condensation = diagnostics.vapor.net_condensation  # cmeout
    λr = incloud.rain.λ  # lamr
    qv_after_step = scratch.number_before  # dum_2D, here the vapor after the step
    T_estimate = scratch.temperature_estimate  # ttmpA
    qsat_adjusted = scratch.q_sat_adjusted  # qvnA
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
    Δt = work.inputs.time_step
    inv_Δt = 1.0 / Δt
    Dcs = settings.ice_to_snow_size
    columns, levels = grid_size(work)

    for k in 1:levels, i in 1:columns
        qc_sed[i, k] = max(qc[i, k] + qc_tendency[i, k] * Δt, 0.0)
        qi_sed[i, k] = max(qi[i, k] + qi_tendency[i, k] * Δt, 0.0)
        nc_sed[i, k] = max(nc[i, k] + nc_tendency[i, k] * Δt, 0.0)
        ni_sed[i, k] = max(ni[i, k] + ni_tendency[i, k] * Δt, 0.0)
        qr_sed[i, k] = max(qr[i, k] + qr_tendency[i, k] * Δt, 0.0)
        nr_sed[i, k] = max(nr[i, k] + nr_tendency[i, k] * Δt, 0.0)
        qs_sed[i, k] = max(qs[i, k] + qs_tendency[i, k] * Δt, 0.0)
        ns_sed[i, k] = max(ns[i, k] + ns_tendency[i, k] * Δt, 0.0)
        qg_sed[i, k] = max(qg[i, k] + qg_tendency[i, k] * Δt, 0.0)
        ng_sed[i, k] = max(ng[i, k] + ng_tendency[i, k] * Δt, 0.0)
        if settings.constant_droplet_number
            nc_sed[i, k] = settings.droplet_number / ρ[i, k] * liquid_fraction[i, k]
        end
        if settings.constant_ice_number
            ni_sed[i, k] = settings.ice_number / ρ[i, k] * ice_fraction[i, k]
        end
        if settings.constant_graupel_number
            ng_sed[i, k] = settings.graupel_number / ρ[i, k] * precip_fraction[i, k]
        end
        if qc_sed[i, k] < smallest_mixing_ratio
            nc_sed[i, k] = 0.0
        end
        if qi_sed[i, k] < smallest_mixing_ratio
            ni_sed[i, k] = 0.0
        end
        if qr_sed[i, k] < smallest_mixing_ratio
            nr_sed[i, k] = 0.0
        end
        if qs_sed[i, k] < smallest_mixing_ratio
            ns_sed[i, k] = 0.0
        end
        if qg_sed[i, k] < smallest_mixing_ratio
            ng_sed[i, k] = 0.0
        end

        # ❄️💧 melting of snow at +2 °C
        if T[i, k] + heating[i, k] / cp * Δt > T_snow_melt
            if qs_sed[i, k] > 0.0
                ΔT_all = -Lf / cp * qs_sed[i, k]  # dum, if all of it changed phase
                if T[i, k] + heating[i, k] / cp * Δt + ΔT_all < T_snow_melt
                    melt_fraction = (T[i, k] + heating[i, k] / cp * Δt - T_snow_melt) * cp / Lf
                    melt_fraction = melt_fraction / qs_sed[i, k]
                    melt_fraction = max(0.0, melt_fraction)
                    melt_fraction = min(1.0, melt_fraction)
                else
                    melt_fraction = 1.0
                end
                qs_tendency[i, k] = qs_tendency[i, k] - melt_fraction * qs_sed[i, k] * inv_Δt
                ns_tendency[i, k] = ns_tendency[i, k] - melt_fraction * ns_sed[i, k] * inv_Δt
                qr_tendency[i, k] = qr_tendency[i, k] + melt_fraction * qs_sed[i, k] * inv_Δt
                nr_tendency[i, k] = nr_tendency[i, k] + melt_fraction * ns_sed[i, k] * inv_Δt
                latent_heating = -Lf * melt_fraction * qs_sed[i, k] * inv_Δt
                heating[i, k] = heating[i, k] + latent_heating
                heating_snow_melting[i, k] = heating_snow_melting[i, k] + latent_heating
            end
        end

        # 🟤💧 melting of graupel at +2 °C
        if T[i, k] + heating[i, k] / cp * Δt > T_snow_melt
            if qg_sed[i, k] > 0.0
                ΔT_all = -Lf / cp * qg_sed[i, k]  # dum, if all of it changed phase
                if T[i, k] + heating[i, k] / cp * Δt + ΔT_all < T_snow_melt
                    melt_fraction = (T[i, k] + heating[i, k] / cp * Δt - T_snow_melt) * cp / Lf
                    melt_fraction = melt_fraction / qg_sed[i, k]
                    melt_fraction = max(0.0, melt_fraction)
                    melt_fraction = min(1.0, melt_fraction)
                else
                    melt_fraction = 1.0
                end
                qg_tendency[i, k] = qg_tendency[i, k] - melt_fraction * qg_sed[i, k] * inv_Δt
                ng_tendency[i, k] = ng_tendency[i, k] - melt_fraction * ng_sed[i, k] * inv_Δt
                qr_tendency[i, k] = qr_tendency[i, k] + melt_fraction * qg_sed[i, k] * inv_Δt
                nr_tendency[i, k] = nr_tendency[i, k] + melt_fraction * ng_sed[i, k] * inv_Δt
                latent_heating = -Lf * melt_fraction * qg_sed[i, k] * inv_Δt
                heating[i, k] = heating[i, k] + latent_heating
                heating_snow_melting[i, k] = heating_snow_melting[i, k] + latent_heating
                if settings.constant_graupel_number && (ng[i, k] + ng_tendency[i, k] * Δt) < 0.0
                    ng_tendency[i, k] = -ng[i, k] * inv_Δt
                end
            end
        end
    end

    # 🌧️ the rain size uses the grid-mean rain (FINDINGS.md)
    for k in 1:levels, i in 1:columns
        rain = size_distribution(hydrometeors.rain, qr_sed[i, k], nr_sed[i, k])
        λr[i, k], nr_sed[i, k] = rain.slope, rain.number
    end

    # 🌧️🧊 freezing of rain at -5 °C
    for k in 1:levels, i in 1:columns
        if T[i, k] + heating[i, k] / cp * Δt < T_rain_freeze
            if qr_sed[i, k] > 0.0
                ΔT_all = Lf / cp * qr_sed[i, k]  # dum, if all of it changed phase
                if T[i, k] + heating[i, k] / cp * Δt + ΔT_all > T_rain_freeze
                    freeze_fraction = -(T[i, k] + heating[i, k] / cp * Δt - T_rain_freeze) * cp / Lf
                    freeze_fraction = freeze_fraction / qr_sed[i, k]
                    freeze_fraction = max(0.0, freeze_fraction)
                    freeze_fraction = min(1.0, freeze_fraction)
                else
                    freeze_fraction = 1.0
                end
                qr_tendency[i, k] = qr_tendency[i, k] - freeze_fraction * qr_sed[i, k] * inv_Δt
                nr_tendency[i, k] = nr_tendency[i, k] - freeze_fraction * nr_sed[i, k] * inv_Δt
                if λr[i, k] < 1.0 / Dcs
                    if settings.do_hail || settings.do_graupel
                        qg_tendency[i, k] = qg_tendency[i, k] + freeze_fraction * qr_sed[i, k] * inv_Δt
                        ng_tendency[i, k] = ng_tendency[i, k] + freeze_fraction * nr_sed[i, k] * inv_Δt
                    else
                        qs_tendency[i, k] = qs_tendency[i, k] + freeze_fraction * qr_sed[i, k] * inv_Δt
                        ns_tendency[i, k] = ns_tendency[i, k] + freeze_fraction * nr_sed[i, k] * inv_Δt
                    end
                else
                    qi_tendency[i, k] = qi_tendency[i, k] + freeze_fraction * qr_sed[i, k] * inv_Δt
                    ni_tendency[i, k] = ni_tendency[i, k] + freeze_fraction * nr_sed[i, k] * inv_Δt
                end
                latent_heating = Lf * freeze_fraction * qr_sed[i, k] * inv_Δt
                heating_rain_freezing[i, k] = heating_rain_freezing[i, k] + latent_heating
                heating[i, k] = heating[i, k] + latent_heating
            end
        end
    end

    settings.prognostic_cloud_ice || return nothing
    homog_size = settings.homogeneous_freezing_size
    max_nicons = settings.maximum_ice_number
    for k in 1:levels, i in 1:columns
        # 🧊💧 melting of cloud ice
        if T[i, k] + heating[i, k] / cp * Δt > T₀
            if qi_sed[i, k] > 0.0
                ΔT_all = -qi_sed[i, k] * Lf / cp  # dum, if all of it changed phase
                if T[i, k] + heating[i, k] / cp * Δt + ΔT_all < T₀
                    melt_fraction = (T[i, k] + heating[i, k] / cp * Δt - T₀) * cp / Lf
                    melt_fraction = melt_fraction / qi_sed[i, k]
                    melt_fraction = max(0.0, melt_fraction)
                    melt_fraction = min(1.0, melt_fraction)
                else
                    melt_fraction = 1.0
                end
                qc_tendency[i, k] = qc_tendency[i, k] + melt_fraction * qi_sed[i, k] * inv_Δt
                qi_melting_mean[i, k] = melt_fraction * qi_sed[i, k] * inv_Δt
                ni_melting_mean[i, k] = 3.0 * melt_fraction * qi_sed[i, k] * inv_Δt /
                                   (4.0 * π * 5.12e-16 * ρ_water)
                nc_tendency[i, k] = nc_tendency[i, k] + ni_melting_mean[i, k]
                qi_tendency[i, k] = ((1.0 - melt_fraction) * qi_sed[i, k] - qi[i, k]) * inv_Δt
                ni_tendency[i, k] = ((1.0 - melt_fraction) * ni_sed[i, k] - ni[i, k]) * inv_Δt
                heating[i, k] = heating[i, k] - Lf * melt_fraction * qi_sed[i, k] * inv_Δt
            end
        end

        # 💧➡️🧊 homogeneous freezing of droplets at -40 °C
        if T[i, k] + heating[i, k] / cp * Δt < 233.15
            if qc_sed[i, k] > 0.0
                ΔT_all = qc_sed[i, k] * Lf / cp  # dum, if all of it changed phase
                if T[i, k] + heating[i, k] / cp * Δt + ΔT_all > 233.15
                    freeze_fraction = -(T[i, k] + heating[i, k] / cp * Δt - 233.15) * cp / Lf
                    freeze_fraction = freeze_fraction / qc_sed[i, k]
                    freeze_fraction = max(0.0, freeze_fraction)
                    freeze_fraction = min(1.0, freeze_fraction)
                else
                    freeze_fraction = 1.0
                end
                qi_tendency[i, k] = qi_tendency[i, k] + freeze_fraction * qc_sed[i, k] * inv_Δt
                qc_homogeneous_freezing_mean[i, k] = freeze_fraction * qc_sed[i, k] * inv_Δt
                nc_homogeneous_freezing_mean[i, k] = freeze_fraction * 3.0 * qc_sed[i, k] /
                                   (4.0 * 3.14 * homog_size^3.0 * 500.0) * inv_Δt
                ni_tendency[i, k] = ni_tendency[i, k] + nc_homogeneous_freezing_mean[i, k]
                qc_tendency[i, k] = ((1.0 - freeze_fraction) * qc_sed[i, k] - qc[i, k]) * inv_Δt
                nc_tendency[i, k] = ((1.0 - freeze_fraction) * nc_sed[i, k] - nc[i, k]) * inv_Δt
                heating[i, k] = heating[i, k] + Lf * freeze_fraction * qc_sed[i, k] * inv_Δt
            end
        end

        # 🔢 ice number limiter
        if ni_tendency[i, k] > 0.0 &&
           ni[i, k] + ni_tendency[i, k] * Δt > max_nicons * ice_fraction[i, k] / ρ[i, k]
            ni_tendency[i, k] = max(0.0, (max_nicons * ice_fraction[i, k] / ρ[i, k] - ni[i, k]) / Δt)
        end

        qv_after_step[i, k] = qv[i, k] + qv_tendency[i, k] * Δt
        T_estimate[i, k] = T[i, k] + heating[i, k] / cp * Δt
    end

    for k in 1:levels, i in 1:columns
        qsat_adjusted[i, k] = saturation_over_water(saturation, T_estimate[i, k], p[i, k]).specific_humidity
    end

    # 🌫️ remove the supersaturation left by the sum of all processes
    for k in 1:levels, i in 1:columns
        if qv_after_step[i, k] > qsat_adjusted[i, k] && qsat_adjusted[i, k] > 0 && settings.remove_supersaturation
            condensation_liquid = (qv_after_step[i, k] - qsat_adjusted[i, k]) /
                  (1.0 + Lv_squared * qsat_adjusted[i, k] / (cp * Rv * (T_estimate[i, k] * T_estimate[i, k]))) *
                  inv_Δt
            net_condensation[i, k] = net_condensation[i, k] + condensation_liquid
            if T_estimate[i, k] > 268.15
                ice_share = 0.0
            elseif T_estimate[i, k] < 238.15
                ice_share = 1.0
            else
                ice_share = (268.15 - T_estimate[i, k]) / 30.0
            end
            latent = Ls * ice_share + Lv * (1.0 - ice_share)
            condensation = (qv_after_step[i, k] - qsat_adjusted[i, k]) /
                  (1.0 + latent * latent * qsat_adjusted[i, k] / (cp * Rv * (T_estimate[i, k] * T_estimate[i, k]))) *
                  inv_Δt
            qc_tendency[i, k] = qc_tendency[i, k] + condensation * (1.0 - ice_share)
            qc_residual_condensation[i, k] = condensation * (1.0 - ice_share)
            qi_tendency[i, k] = qi_tendency[i, k] + condensation * ice_share
            qi_residual_deposition[i, k] = condensation * ice_share
            qv_tendency[i, k] = qv_tendency[i, k] - condensation
            qv_residual_condensation[i, k] = -condensation
            heating[i, k] = heating[i, k] + condensation * (1.0 - ice_share) * Lv + condensation * ice_share * Ls
        end
    end
    return nothing
end
