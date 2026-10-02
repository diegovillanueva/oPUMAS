# ⬇️🌧️ Stage 12, sedimentation (micro_pumas_v1.F90 lines 3507–3730, with the
#    `Sedimentation`, `Sedimentation_implicit` and `implicit_fall` subroutines). The
#    implicit scheme follows Guo et al. (2021); the explicit one uses sub-steps.

"""
    implicit_fall!(q, flux, deltat, i, zint, speed, pdel, levels) -> precipitation

⬇️ Twin of `implicit_fall` for column `i`: updates the mixing ratio `q[i, :]` in place,
writes the cumulative fall-out `flux[i, :]` and returns the surface value.
"""
function implicit_fall!(q, flux, deltat, i, zint, speed, pdel, levels)
    dz = Vector{Float64}(undef, levels)
    dd = Vector{Float64}(undef, levels)
    qm = Vector{Float64}(undef, levels)
    for k in 1:levels
        dz[k] = zint[i, k] - zint[i, k+1]
        dd[k] = deltat * speed[i, k]
        q[i, k] = q[i, k] * pdel[i, k]
    end
    qm[1] = q[i, 1] / (dz[1] + dd[1])
    for k in 2:levels
        qm[k] = (q[i, k] + dd[k-1] * qm[k-1]) / (dz[k] + dd[k])
    end
    for k in 1:levels
        qm[k] = qm[k] * dz[k]
    end
    flux[i, 1] = q[i, 1] - qm[1]
    for k in 2:levels
        flux[i, k] = flux[i, k-1] + q[i, k] - qm[k]
    end
    precipitation = flux[i, levels]
    for k in 1:levels
        q[i, k] = qm[k] / pdel[i, k]
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
    (; g) = host_short_names(parameters)
    F = work.fortran
    (; zint, pdel) = F
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)
    after = copy(amount)
    falling = similar(amount)
    precipitation = Vector{Float64}(undef, columns)
    for i in 1:columns
        precipitation[i] = implicit_fall!(after, falling, deltat, i, zint, speed, pdel,
            levels)
    end
    for k in 1:levels, i in 1:columns
        if check_qsmall
            if falling[i, k] >= qsmall && flux !== nothing
                flux[i, k+1] = flux[i, k+1] + falling[i, k] / g / deltat
            end
        elseif flux !== nothing
            flux[i, k+1] = flux[i, k+1] + falling[i, k] / g / deltat
        end
        if sedimentation !== nothing
            sedimentation[i, k] = sedimentation[i, k] + (after[i, k] - amount[i, k]) / deltat
        end
        tendency[i, k] = tendency[i, k] + (after[i, k] - amount[i, k]) / deltat
    end
    for i in 1:columns
        if precipitation[i] >= 0.0
            if surface !== nothing
                surface[i] = surface[i] + precipitation[i] / g / deltat / 1000.0
            end
            if surface_ice !== nothing
                surface_ice[i] = surface_ice[i] + precipitation[i] / g / deltat / 1000.0
            end
        end
    end
    return nothing
end

"""
    sub_steps(mass_speed, number_speed, work)

⏱️ Explicit sub-steps per column, `nstep` and `rnstep`, so that no level loses more than
its content in one sub-step.
"""
function sub_steps(mass_speed, number_speed, work)
    pdel_inv = work.fortran.pdel_inv
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)
    nstep = Vector{Int}(undef, columns)
    rnstep = Vector{Float64}(undef, columns)
    for i in 1:columns
        mass_max = maximum(mass_speed[i, k] * pdel_inv[i, k] for k in 1:levels)
        number_max = maximum(number_speed[i, k] * pdel_inv[i, k] for k in 1:levels)
        nstep[i] = 1 + trunc(Int, max(mass_max, number_max) * deltat)
        rnstep[i] = 1.0 / Float64(Float32(nstep[i]))
    end
    return nstep, rnstep
end

"""
    sediment_explicit!(tendency, amount, speed, nstep, rnstep, work, parameters,
                       falls; sedimentation, surface, flux, latent_heat, evaporation,
                       heating, vapor, fraction, surface_ice)

⬇️ Twin of `Sedimentation`: `nstep` forward sub-steps of the flux divergence. With a
cloud `fraction`, the flux entering the clear part of a level evaporates at once.
"""
function sediment_explicit!(tendency, amount, speed, nstep, rnstep, work, parameters,
    falls; sedimentation = nothing, surface = nothing, flux = nothing,
    latent_heat = nothing, evaporation = nothing, heating = nothing, vapor = nothing,
    fraction = nothing, surface_ice = nothing)
    (; g) = host_short_names(parameters)
    pdel_inv = work.fortran.pdel_inv
    deltat = work.inputs.time_step
    columns, levels = grid_size(work)
    dum1 = Vector{Float64}(undef, levels)
    fallout = Vector{Float64}(undef, levels + 1)   # fallout[k + 1] is faloutx(i, k)
    for i in 1:columns
        nstepmax = nstep[i]
        rnstepmax = rnstep[i]
        dum1[1] = 0.0
        if fraction !== nothing
            for k in 2:levels
                dum1[k] = fraction[i, k] / fraction[i, k-1]
                dum1[k] = min(dum1[k], 1.0)
            end
        else
            for k in 2:levels
                dum1[k] = 1.0
            end
        end
        for n in 1:nstepmax
            fallout[1] = 0.0
            for k in 1:levels
                fallout[k+1] = falls ? speed[i, k] * amount[i, k] : 0.0
            end
            for k in 1:levels
                faltndx = (fallout[k+1] - dum1[k] * fallout[k]) * pdel_inv[i, k]
                tendency[i, k] = tendency[i, k] - faltndx * rnstepmax
                if sedimentation !== nothing
                    sedimentation[i, k] = sedimentation[i, k] - faltndx * rnstepmax
                end
                amount[i, k] = amount[i, k] - faltndx * deltat * rnstepmax
                if k > 1 && (evaporation !== nothing || vapor !== nothing ||
                             heating !== nothing)
                    faltndqxe = (fallout[k+1] - fallout[k]) * pdel_inv[i, k]
                    if evaporation !== nothing
                        evaporation[i, k] = evaporation[i, k] -
                                            (faltndqxe - faltndx) * rnstepmax
                    end
                    if vapor !== nothing
                        vapor[i, k] = vapor[i, k] - (faltndqxe - faltndx) * rnstepmax
                    end
                    if heating !== nothing
                        heating[i, k] = heating[i, k] +
                                        (faltndqxe - faltndx) * latent_heat * rnstepmax
                    end
                end
                if flux !== nothing
                    flux[i, k+1] = flux[i, k+1] + fallout[k+1] / g * rnstepmax
                end
            end
            if surface !== nothing
                surface[i] = surface[i] + fallout[levels+1] / g * rnstepmax / 1000.0
            end
            if surface_ice !== nothing
                surface_ice[i] = surface_ice[i] + fallout[levels+1] / g * rnstepmax / 1000.0
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
    (; g, xxlv, xxls) = host_short_names(parameters)
    F = work.fortran
    (; rho, lcldm, icldm, tlat, qvlat, prect, preci) = F
    (; qctend, nctend, qitend, nitend, qrtend, nrtend, qstend, nstend, qgtend, ngtend) = F
    (; dumc, dumnc, dumi, dumni, dumr, dumnr, dums, dumns, dumg, dumng) = F
    (; fc, fnc, fi, fni, fr, fnr, fs, fns, fg, fng) = F
    (; lflx, iflx, rflx, sflx, gflx, tlat_i, tlat_l, qvlat_i, qvlat_l) = F
    (; prect_l, prect_i, prect_r, prect_s, prect_g, preci_i, preci_s, preci_g) = F
    columns, levels = grid_size(work)

    if settings.implicit_sedimentation
        vfactor, vfac_drop, vfac_ice = 1.0, 1.0, 1.0
        for (factor, speeds) in ((vfac_drop, (fc, fnc)), (vfac_ice, (fi, fni)),
            (vfactor, (fr, fnr)), (vfactor, (fs, fns)), (vfactor, (fg, fng)))
            for speed in speeds, k in 1:levels, i in 1:columns
                speed[i, k] = factor * speed[i, k] / g / rho[i, k]
            end
        end
        sediment_implicit!(qctend, dumc, fc, work, parameters, false; flux = lflx,
            sedimentation = F.qcsedten, surface = prect_l)
        sediment_implicit!(nctend, dumnc, fnc, work, parameters, false;
            sedimentation = F.ncsedten)
        sediment_implicit!(qitend, dumi, fi, work, parameters, false; flux = iflx,
            sedimentation = F.qisedten, surface = prect_i, surface_ice = preci_i)
        sediment_implicit!(nitend, dumni, fni, work, parameters, false;
            sedimentation = F.nisedten)
        sediment_implicit!(qrtend, dumr, fr, work, parameters, true; flux = rflx,
            sedimentation = F.qrsedten, surface = prect_r)
        sediment_implicit!(nrtend, dumnr, fnr, work, parameters, true;
            sedimentation = F.nrsedten)
        sediment_implicit!(qstend, dums, fs, work, parameters, true; flux = sflx,
            sedimentation = F.qssedten, surface = prect_s, surface_ice = preci_s)
        sediment_implicit!(nstend, dumns, fns, work, parameters, true;
            sedimentation = F.nssedten)
        sediment_implicit!(qgtend, dumg, fg, work, parameters, true; flux = gflx,
            sedimentation = F.qgsedten, surface = prect_g, surface_ice = preci_g)
        sediment_implicit!(ngtend, dumng, fng, work, parameters, true;
            sedimentation = F.ngsedten)
    else
        ice_on = settings.prognostic_cloud_ice
        nstep, rnstep = sub_steps(fi, fni, work)
        sediment_explicit!(qitend, dumi, fi, nstep, rnstep, work, parameters, ice_on;
            sedimentation = F.qisedten, surface = prect_i, flux = iflx,
            latent_heat = xxls, evaporation = F.qisevap, heating = tlat_i,
            vapor = qvlat_i, fraction = icldm, surface_ice = preci_i)
        sediment_explicit!(nitend, dumni, fni, nstep, rnstep, work, parameters, ice_on;
            fraction = icldm, sedimentation = F.nisedten)
        nstep, rnstep = sub_steps(fc, fnc, work)
        sediment_explicit!(qctend, dumc, fc, nstep, rnstep, work, parameters, true;
            sedimentation = F.qcsedten, surface = prect_l, flux = lflx,
            latent_heat = xxlv, evaporation = F.qcsevap, heating = tlat_l,
            vapor = qvlat_l, fraction = lcldm)
        sediment_explicit!(nctend, dumnc, fnc, nstep, rnstep, work, parameters, true;
            fraction = lcldm, sedimentation = F.ncsedten)
        nstep, rnstep = sub_steps(fr, fnr, work)
        sediment_explicit!(qrtend, dumr, fr, nstep, rnstep, work, parameters, true;
            sedimentation = F.qrsedten, surface = prect_r, flux = rflx)
        sediment_explicit!(nrtend, dumnr, fnr, nstep, rnstep, work, parameters, true;
            sedimentation = F.nrsedten)
        nstep, rnstep = sub_steps(fs, fns, work)
        sediment_explicit!(qstend, dums, fs, nstep, rnstep, work, parameters, true;
            sedimentation = F.qssedten, surface = prect_s, flux = sflx,
            surface_ice = preci_s)
        sediment_explicit!(nstend, dumns, fns, nstep, rnstep, work, parameters, true;
            sedimentation = F.nssedten)
        nstep, rnstep = sub_steps(fg, fng, work)
        sediment_explicit!(qgtend, dumg, fg, nstep, rnstep, work, parameters, true;
            sedimentation = F.qgsedten, surface = prect_g, flux = gflx,
            surface_ice = preci_g)
        sediment_explicit!(ngtend, dumng, fng, nstep, rnstep, work, parameters, true;
            sedimentation = F.ngsedten)
    end

    for k in 1:levels, i in 1:columns
        tlat[i, k] = tlat[i, k] + tlat_i[i, k] + tlat_l[i, k]
        qvlat[i, k] = qvlat[i, k] + qvlat_i[i, k] + qvlat_l[i, k]
    end
    for i in 1:columns
        prect[i] = prect[i] + prect_i[i] + prect_l[i] + prect_r[i] + prect_s[i] +
                   prect_g[i]
        preci[i] = preci[i] + preci_i[i] + preci_s[i] + preci_g[i]
    end
    return nothing
end
