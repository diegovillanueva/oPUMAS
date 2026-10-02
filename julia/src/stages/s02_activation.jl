# 🌱❄️ Stage 2, activation and ice nucleation (micro_pumas_v1.F90 lines 1715–1783).

"""
    activation_and_nucleation!(work, parameters)

🌱 Stage 2. Droplet activation from the host's `npccn` is added to `nc` at once; ice
nucleates from the activated nuclei `naai` below -5 °C (with RH above 105 % over ice
unless `ice_nucleation_without_rh`).
"""
function activation_and_nucleation!(work, parameters)
    (; settings) = parameters
    (; T_ice_nucleation) = host_constants(parameters)
    (; inputs, state, incloud, rates, diagnostics) = work
    T = inputs.T  # t
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    droplet_activation = inputs.activated_droplets  # npccn
    ice_nuclei = inputs.activated_ice_nuclei  # naai
    ρ = state.air.density  # rho
    ice_fraction = incloud.ice_fraction  # icldm
    rh = state.air.relative_humidity  # relhum
    e_sat_liquid = state.air.vapor_pressure_liquid  # esl
    e_sat_ice = state.air.vapor_pressure_ice  # esi
    nc_activated = diagnostics.cloud_liquid.activation  # ncal
    ni_available = diagnostics.cloud_ice.available_nuclei  # ncai
    ni_nucleation = rates.cloud_ice.nucleation_from_vapor.number  # nnuccd
    qi_nucleation = rates.cloud_ice.nucleation_from_vapor.mass  # mnuccd
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= smallest_mixing_ratio
            nc[i, k] = max(nc[i, k] + droplet_activation[i, k] * Δt, 0.0)
            nc_activated[i, k] = droplet_activation[i, k]
        else
            nc_activated[i, k] = 0.0
        end
        if T[i, k] < T_ice_nucleation
            ni_available[i, k] = ice_nuclei[i, k] * Δt * ρ[i, k]
        else
            ni_available[i, k] = 0.0
        end
        if settings.prognostic_cloud_ice
            nucleates = ice_nuclei[i, k] > 0.0 && T[i, k] < T_ice_nucleation
            if !settings.ice_nucleation_without_rh
                nucleates = nucleates &&
                            rh[i, k] * e_sat_liquid[i, k] / e_sat_ice[i, k] > 1.05
            end
            if nucleates
                ni_nucleation[i, k] = ice_nuclei[i, k] * ice_fraction[i, k]
                ni_nucleation[i, k] = max(ni_nucleation[i, k], 0.0)
                qi_nucleation[i, k] = ni_nucleation[i, k] * new_ice_crystal_mass
            else
                ni_nucleation[i, k] = 0.0
                qi_nucleation[i, k] = 0.0
            end
        end
    end
    return nothing
end
