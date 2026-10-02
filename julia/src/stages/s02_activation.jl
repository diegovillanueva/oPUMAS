# 🌱❄️ Stage 2, activation and ice nucleation (micro_pumas_v1.F90 lines 1715–1783).

"""
    activation_and_nucleation!(work, parameters)

🌱 Stage 2. Droplet activation from the host's `npccn` is added to `nc` at once; ice
nucleates from the activated nuclei `naai` below -5 °C (with RH above 105 % over ice
unless `ice_nucleation_without_rh`).
"""
function activation_and_nucleation!(work, parameters)
    (; settings) = parameters
    (; icenuct) = host_short_names(parameters)
    F = work.fortran
    (; t, qc, nc, npccn, naai, rho, icldm, relhum, esl, esi, ncal, ncai) = F
    (; nnuccd, mnuccd) = F
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= qsmall
            nc[i, k] = max(nc[i, k] + npccn[i, k] * deltat, 0.0)
            ncal[i, k] = npccn[i, k]
        else
            ncal[i, k] = 0.0
        end
        if t[i, k] < icenuct
            ncai[i, k] = naai[i, k] * deltat * rho[i, k]
        else
            ncai[i, k] = 0.0
        end
        if settings.prognostic_cloud_ice
            nucleates = naai[i, k] > 0.0 && t[i, k] < icenuct
            if !settings.ice_nucleation_without_rh
                nucleates = nucleates && relhum[i, k] * esl[i, k] / esi[i, k] > 1.05
            end
            if nucleates
                nnuccd[i, k] = naai[i, k] * icldm[i, k]
                nnuccd[i, k] = max(nnuccd[i, k], 0.0)
                mnuccd[i, k] = nnuccd[i, k] * mi0
            else
                nnuccd[i, k] = 0.0
                mnuccd[i, k] = 0.0
            end
        end
    end
    return nothing
end
