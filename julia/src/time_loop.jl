# ⏱️📈 Offline single-column time loop: repeated steps with the host forcing held fixed
#    (twin of the `nsteps` loop of tests/column/driver.F90).

"""
    advance_state!(work, parameters)

⏱️ Forward Euler with the tendencies of the last step, as tests/column/driver.F90 does
between two calls: `T += heating / cpair * Δt`, and `x += tendency * Δt` for vapor and
every q and N, with negative values set to zero. This is the only place that changes
`work.inputs`; everything else the host provides (cloud fractions, nuclei, pressure)
stays as it was read.

Fortran twin: `advance_state` of tests/column/driver.F90.
"""
function advance_state!(work, parameters)
    (; inputs, tendencies) = work
    (; cp) = host_constants(parameters)
    Δt = inputs.time_step
    @. inputs.T = inputs.T + tendencies.heating / cp * Δt
    forward_euler!(inputs.q_vapor, tendencies.vapor, Δt)
    for species in (:cloud_liquid, :cloud_ice, :rain, :snow, :graupel), moment in (:q, :N)
        forward_euler!(inputs[species][moment], tendencies[species][moment], Δt)
    end
    return nothing
end

"🧮 `field = field + tendency * Δt`, then zero where negative (NaN stays, as in `where`)."
function forward_euler!(field, tendency, Δt)
    for i in eachindex(field, tendency)
        field[i] = field[i] + tendency[i] * Δt
        field[i] < 0.0 && (field[i] = 0.0)
    end
    return field
end

"""
    time_loop!(work, parameters, steps; after_step = (work, step) -> nothing) -> work

⏱️ `steps` calls of [`microphysics_step!`](@ref), with [`advance_state!`](@ref) between
two calls. `after_step` sees the workspace after each step, e.g. to record profiles.

Fortran twin: the `nsteps` loop of `pumas_column_driver` (tests/column/driver.F90).
"""
function time_loop!(work, parameters, steps; after_step = (work, step) -> nothing)
    for step in 1:steps
        step > 1 && advance_state!(work, parameters)
        microphysics_step!(work, parameters)
        after_step(work, step)
    end
    return work
end
