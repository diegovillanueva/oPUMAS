# ⬇️🌧️ Stage 12, sedimentation (micro_pumas_v1.F90 lines 3507–3730, with the
#    `Sedimentation`, `Sedimentation_implicit` and `implicit_fall` subroutines). The
#    implicit scheme follows Guo et al. (2021); the explicit one uses sub-steps.

"""
    implicit_fall!(q, flux, deltat, i, zint, speed, pdel, levels) -> precipitation

⬇️ Twin of `implicit_fall` for column `i`: updates the mixing ratio `q[i, :]` in place,
writes the cumulative fall-out `flux[i, :]` and returns the surface value.
"""
function implicit_fall!(q, flux, Δt, i, zint, speed, pdel, levels)
    layer_depth = Vector{Float64}(undef, levels)
    fall_depth = Vector{Float64}(undef, levels)
    q_new = Vector{Float64}(undef, levels)
    for k in 1:levels
        layer_depth[k] = zint[i, k] - zint[i, k + 1]
        fall_depth[k] = Δt * speed[i, k]
        q[i, k] = q[i, k] * pdel[i, k]
    end
    q_new[1] = q[i, 1] / (layer_depth[1] + fall_depth[1])
    for k in 2:levels
        q_new[k] = (q[i, k] + fall_depth[k - 1] * q_new[k - 1]) /
                   (layer_depth[k] + fall_depth[k])
    end
    for k in 1:levels
        q_new[k] = q_new[k] * layer_depth[k]
    end
    flux[i, 1] = q[i, 1] - q_new[1]
    for k in 2:levels
        flux[i, k] = flux[i, k - 1] + q[i, k] - q_new[k]
    end
    precipitation = flux[i, levels]
    for k in 1:levels
        q[i, k] = q_new[k] / pdel[i, k]
    end
    return precipitation
end

"""
    sediment_implicit!(tendency, amount, speed, work, parameters, check_qsmall;
                       flux, sedimentation, surface, surface_ice)

⬇️ Twin of `Sedimentation_implicit`. Optional outputs are `nothing` when the Fortran
call leaves them out.
"""
function sediment_implicit!(tendency, amount, speed, work, parameters, check_qsmall;
        flux = nothing, sedimentation = nothing, surface = nothing, surface_ice = nothing)
    (; g) = host_constants(parameters)
    (; inputs, scratch) = work
    z_interface = scratch.interface_height  # zint
    Δp = inputs.pressure_thickness  # pdel
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)
    after = copy(amount)
    falling = similar(amount)
    precipitation = Vector{Float64}(undef, columns)
    for i in 1:columns
        precipitation[i] = implicit_fall!(after, falling, Δt, i, z_interface, speed, Δp,
            levels)
    end
    for k in 1:levels, i in 1:columns
        if check_qsmall
            if falling[i, k] >= smallest_mixing_ratio && flux !== nothing
                flux[i, k + 1] = flux[i, k + 1] + falling[i, k] / g / Δt
            end
        elseif flux !== nothing
            flux[i, k + 1] = flux[i, k + 1] + falling[i, k] / g / Δt
        end
        if sedimentation !== nothing
            sedimentation[i, k] = sedimentation[i, k] + (after[i, k] - amount[i, k]) / Δt
        end
        tendency[i, k] = tendency[i, k] + (after[i, k] - amount[i, k]) / Δt
    end
    for i in 1:columns
        if precipitation[i] >= 0.0
            if surface !== nothing
                surface[i] = surface[i] + precipitation[i] / g / Δt / 1000.0
            end
            if surface_ice !== nothing
                surface_ice[i] = surface_ice[i] + precipitation[i] / g / Δt / 1000.0
            end
        end
    end
    return nothing
end

"""
    sub_steps(mass_speed, number_speed, work)

⏱️ Explicit sub-steps per column, `substeps` and `inv_substeps` (`nstep`, `rnstep`), so
that no level loses more than its content in one sub-step.
"""
function sub_steps(mass_speed, number_speed, work)
    (; scratch) = work
    inv_Δp = scratch.inverse_pressure_thickness  # pdel_inv
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)
    substeps = Vector{Int}(undef, columns)
    inv_substeps = Vector{Float64}(undef, columns)
    for i in 1:columns
        mass_max = maximum(mass_speed[i, k] * inv_Δp[i, k] for k in 1:levels)
        number_max = maximum(number_speed[i, k] * inv_Δp[i, k] for k in 1:levels)
        substeps[i] = 1 + trunc(Int, max(mass_max, number_max) * Δt)
        inv_substeps[i] = 1.0 / Float64(Float32(substeps[i]))
    end
    return substeps, inv_substeps
end

"""
    sediment_explicit!(tendency, amount, speed, substeps, inv_substeps, work, parameters,
                       falls; sedimentation, surface, flux, latent_heat, evaporation,
                       heating, vapor, fraction, surface_ice)

⬇️ Twin of `Sedimentation`: `substeps` forward sub-steps of the flux divergence. With a
cloud `fraction`, the flux entering the clear part of a level evaporates at once.
"""
function sediment_explicit!(
        tendency, amount, speed, substeps, inv_substeps, work, parameters,
        falls; sedimentation = nothing, surface = nothing, flux = nothing,
        latent_heat = nothing, evaporation = nothing, heating = nothing, vapor = nothing,
        fraction = nothing, surface_ice = nothing)
    (; g) = host_constants(parameters)
    (; scratch) = work
    inv_Δp = scratch.inverse_pressure_thickness  # pdel_inv
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)
    inflow_share = Vector{Float64}(undef, levels)
    fallout = Vector{Float64}(undef, levels + 1)   # fallout[k + 1] is faloutx(i, k)
    for i in 1:columns
        steps = substeps[i]
        inv_steps = inv_substeps[i]
        inflow_share[1] = 0.0
        if fraction !== nothing
            for k in 2:levels
                inflow_share[k] = fraction[i, k] / fraction[i, k - 1]
                inflow_share[k] = min(inflow_share[k], 1.0)
            end
        else
            for k in 2:levels
                inflow_share[k] = 1.0
            end
        end
        for n in 1:steps
            fallout[1] = 0.0
            for k in 1:levels
                fallout[k + 1] = falls ? speed[i, k] * amount[i, k] : 0.0
            end
            for k in 1:levels
                flux_divergence = (fallout[k + 1] - inflow_share[k] * fallout[k]) *
                                  inv_Δp[i, k]
                tendency[i, k] = tendency[i, k] - flux_divergence * inv_steps
                if sedimentation !== nothing
                    sedimentation[i, k] = sedimentation[i, k] - flux_divergence * inv_steps
                end
                amount[i, k] = amount[i, k] - flux_divergence * Δt * inv_steps
                if k > 1 && (evaporation !== nothing || vapor !== nothing ||
                    heating !== nothing)
                    flux_divergence_unscaled = (fallout[k + 1] - fallout[k]) * inv_Δp[i, k]
                    if evaporation !== nothing
                        evaporation[i, k] = evaporation[i, k] -
                                            (flux_divergence_unscaled - flux_divergence) *
                                            inv_steps
                    end
                    if vapor !== nothing
                        vapor[i, k] = vapor[i, k] -
                                      (flux_divergence_unscaled - flux_divergence) *
                                      inv_steps
                    end
                    if heating !== nothing
                        heating[i, k] = heating[i, k] +
                                        (flux_divergence_unscaled - flux_divergence) *
                                        latent_heat * inv_steps
                    end
                end
                if flux !== nothing
                    flux[i, k + 1] = flux[i, k + 1] + fallout[k + 1] / g * inv_steps
                end
            end
            if surface !== nothing
                surface[i] = surface[i] + fallout[levels + 1] / g * inv_steps / 1000.0
            end
            if surface_ice !== nothing
                surface_ice[i] = surface_ice[i] +
                                 fallout[levels + 1] / g * inv_steps / 1000.0
            end
        end
    end
    return nothing
end

"""
    sedimentation!(work, parameters)

⬇️ Stage 12. Sediments the mass and number of every species, then sums the cloud
evaporation terms into `tlat`/`qvlat` and the surface rates into `prect`/`preci`.
"""
function sedimentation!(work, parameters)
    (; settings) = parameters
    (; g, Lv, Ls) = host_constants(parameters)
    (; state, incloud, process_rates, tendencies, diagnostics, scratch) = work
    liquid_means = process_rates.cloud_liquid
    ice_means = process_rates.cloud_ice
    rain_means = process_rates.rain
    snow_means = process_rates.snow
    graupel_means = process_rates.graupel
    qc_sedimentation = liquid_means.sedimentation.mass  # qcsedten
    nc_sedimentation = liquid_means.sedimentation.number  # ncsedten
    qi_sedimentation = ice_means.sedimentation.mass  # qisedten
    ni_sedimentation = ice_means.sedimentation.number  # nisedten
    qr_sedimentation = rain_means.sedimentation.mass  # qrsedten
    nr_sedimentation = rain_means.sedimentation.number  # nrsedten
    qs_sedimentation = snow_means.sedimentation.mass  # qssedten
    ns_sedimentation = snow_means.sedimentation.number  # nssedten
    qg_sedimentation = graupel_means.sedimentation.mass  # qgsedten
    ng_sedimentation = graupel_means.sedimentation.number  # ngsedten
    qi_sedimentation_sublimation = ice_means.sedimentation_sublimation.mass  # qisevap
    qc_sedimentation_evaporation = liquid_means.sedimentation_evaporation.mass  # qcsevap
    ρ = state.air.density  # rho
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    heating = tendencies.heating  # tlat
    qv_tendency = tendencies.vapor  # qvlat
    precip_surface = diagnostics.surface.total  # prect
    precip_surface_ice = diagnostics.surface.ice  # preci
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
    fall_qc = scratch.pressure_fall_speed.cloud_liquid.q  # fc
    fall_nc = scratch.pressure_fall_speed.cloud_liquid.N  # fnc
    fall_qi = scratch.pressure_fall_speed.cloud_ice.q  # fi
    fall_ni = scratch.pressure_fall_speed.cloud_ice.N  # fni
    fall_qr = scratch.pressure_fall_speed.rain.q  # fr
    fall_nr = scratch.pressure_fall_speed.rain.N  # fnr
    fall_qs = scratch.pressure_fall_speed.snow.q  # fs
    fall_ns = scratch.pressure_fall_speed.snow.N  # fns
    fall_qg = scratch.pressure_fall_speed.graupel.q  # fg
    fall_ng = scratch.pressure_fall_speed.graupel.N  # fng
    flux_cloud_liquid = diagnostics.flux.cloud_liquid  # lflx
    flux_cloud_ice = diagnostics.flux.cloud_ice  # iflx
    flux_rain = diagnostics.flux.rain  # rflx
    flux_snow = diagnostics.flux.snow  # sflx
    flux_graupel = diagnostics.flux.graupel  # gflx
    heating_sed_ice = scratch.sedimentation_heating.cloud_ice  # tlat_i
    heating_sed_liquid = scratch.sedimentation_heating.cloud_liquid  # tlat_l
    qv_sed_ice = scratch.sedimentation_vapor.cloud_ice  # qvlat_i
    qv_sed_liquid = scratch.sedimentation_vapor.cloud_liquid  # qvlat_l
    surface_cloud_liquid = scratch.surface_precipitation.cloud_liquid  # prect_l
    surface_cloud_ice = scratch.surface_precipitation.cloud_ice  # prect_i
    surface_rain = scratch.surface_precipitation.rain  # prect_r
    surface_snow = scratch.surface_precipitation.snow  # prect_s
    surface_graupel = scratch.surface_precipitation.graupel  # prect_g
    surface_frozen_cloud_ice = scratch.surface_ice.cloud_ice  # preci_i
    surface_frozen_snow = scratch.surface_ice.snow  # preci_s
    surface_frozen_graupel = scratch.surface_ice.graupel  # preci_g
    columns, levels = grid_size(work)

    if settings.implicit_sedimentation
        vfactor, vfac_drop, vfac_ice = 1.0, 1.0, 1.0
        for (factor, speeds) in (
            (vfac_drop, (fall_qc, fall_nc)), (vfac_ice, (fall_qi, fall_ni)),
            (vfactor, (fall_qr, fall_nr)), (vfactor, (fall_qs, fall_ns)), (
                vfactor, (fall_qg, fall_ng)))
            for speed in speeds, k in 1:levels, i in 1:columns
                speed[i, k] = factor * speed[i, k] / g / ρ[i, k]
            end
        end
        sediment_implicit!(
            qc_tendency, qc_sed, fall_qc, work, parameters, false; flux = flux_cloud_liquid,
            sedimentation = qc_sedimentation, surface = surface_cloud_liquid)
        sediment_implicit!(nc_tendency, nc_sed, fall_nc, work, parameters, false;
            sedimentation = nc_sedimentation)
        sediment_implicit!(
            qi_tendency, qi_sed, fall_qi, work, parameters, false; flux = flux_cloud_ice,
            sedimentation = qi_sedimentation, surface = surface_cloud_ice,
            surface_ice = surface_frozen_cloud_ice)
        sediment_implicit!(ni_tendency, ni_sed, fall_ni, work, parameters, false;
            sedimentation = ni_sedimentation)
        sediment_implicit!(
            qr_tendency, qr_sed, fall_qr, work, parameters, true; flux = flux_rain,
            sedimentation = qr_sedimentation, surface = surface_rain)
        sediment_implicit!(nr_tendency, nr_sed, fall_nr, work, parameters, true;
            sedimentation = nr_sedimentation)
        sediment_implicit!(
            qs_tendency, qs_sed, fall_qs, work, parameters, true; flux = flux_snow,
            sedimentation = qs_sedimentation, surface = surface_snow,
            surface_ice = surface_frozen_snow)
        sediment_implicit!(ns_tendency, ns_sed, fall_ns, work, parameters, true;
            sedimentation = ns_sedimentation)
        sediment_implicit!(
            qg_tendency, qg_sed, fall_qg, work, parameters, true; flux = flux_graupel,
            sedimentation = qg_sedimentation, surface = surface_graupel,
            surface_ice = surface_frozen_graupel)
        sediment_implicit!(ng_tendency, ng_sed, fall_ng, work, parameters, true;
            sedimentation = ng_sedimentation)
    else
        ice_on = settings.prognostic_cloud_ice
        substeps, inv_substeps = sub_steps(fall_qi, fall_ni, work)
        sediment_explicit!(
            qi_tendency, qi_sed, fall_qi, substeps, inv_substeps, work, parameters, ice_on;
            sedimentation = qi_sedimentation, surface = surface_cloud_ice,
            flux = flux_cloud_ice,
            latent_heat = Ls, evaporation = qi_sedimentation_sublimation,
            heating = heating_sed_ice,
            vapor = qv_sed_ice, fraction = ice_fraction,
            surface_ice = surface_frozen_cloud_ice)
        sediment_explicit!(
            ni_tendency, ni_sed, fall_ni, substeps, inv_substeps, work, parameters, ice_on;
            fraction = ice_fraction, sedimentation = ni_sedimentation)
        substeps, inv_substeps = sub_steps(fall_qc, fall_nc, work)
        sediment_explicit!(
            qc_tendency, qc_sed, fall_qc, substeps, inv_substeps, work, parameters, true;
            sedimentation = qc_sedimentation, surface = surface_cloud_liquid,
            flux = flux_cloud_liquid,
            latent_heat = Lv, evaporation = qc_sedimentation_evaporation,
            heating = heating_sed_liquid,
            vapor = qv_sed_liquid, fraction = liquid_fraction)
        sediment_explicit!(
            nc_tendency, nc_sed, fall_nc, substeps, inv_substeps, work, parameters, true;
            fraction = liquid_fraction, sedimentation = nc_sedimentation)
        substeps, inv_substeps = sub_steps(fall_qr, fall_nr, work)
        sediment_explicit!(
            qr_tendency, qr_sed, fall_qr, substeps, inv_substeps, work, parameters, true;
            sedimentation = qr_sedimentation, surface = surface_rain, flux = flux_rain)
        sediment_explicit!(
            nr_tendency, nr_sed, fall_nr, substeps, inv_substeps, work, parameters, true;
            sedimentation = nr_sedimentation)
        substeps, inv_substeps = sub_steps(fall_qs, fall_ns, work)
        sediment_explicit!(
            qs_tendency, qs_sed, fall_qs, substeps, inv_substeps, work, parameters, true;
            sedimentation = qs_sedimentation, surface = surface_snow, flux = flux_snow,
            surface_ice = surface_frozen_snow)
        sediment_explicit!(
            ns_tendency, ns_sed, fall_ns, substeps, inv_substeps, work, parameters, true;
            sedimentation = ns_sedimentation)
        substeps, inv_substeps = sub_steps(fall_qg, fall_ng, work)
        sediment_explicit!(
            qg_tendency, qg_sed, fall_qg, substeps, inv_substeps, work, parameters, true;
            sedimentation = qg_sedimentation, surface = surface_graupel,
            flux = flux_graupel,
            surface_ice = surface_frozen_graupel)
        sediment_explicit!(
            ng_tendency, ng_sed, fall_ng, substeps, inv_substeps, work, parameters, true;
            sedimentation = ng_sedimentation)
    end

    for k in 1:levels, i in 1:columns
        heating[i, k] = heating[i, k] + heating_sed_ice[i, k] + heating_sed_liquid[i, k]
        qv_tendency[i, k] = qv_tendency[i, k] + qv_sed_ice[i, k] + qv_sed_liquid[i, k]
    end
    for i in 1:columns
        precip_surface[i] = precip_surface[i] + surface_cloud_ice[i] +
                            surface_cloud_liquid[i] + surface_rain[i] + surface_snow[i] +
                            surface_graupel[i]
        precip_surface_ice[i] = precip_surface_ice[i] + surface_frozen_cloud_ice[i] +
                                surface_frozen_snow[i] + surface_frozen_graupel[i]
    end
    return nothing
end
