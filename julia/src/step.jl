# 🔁🌧️ One time step of the microphysics: the fifteen stages of micro_pumas_tend in
#    Fortran order (julia/DESIGN.md, "Stages").

"""
    microphysics_step!(work, parameters) -> work

🔁 Twin of `micro_pumas_tend`. Reads `work.inputs`, overwrites every other container of
`work` (from [`allocate`](@ref)) with the tendencies, diagnostics and process rates of
one step.
"""
function microphysics_step!(work, parameters)
    setup_step!(work, parameters)                       # 1  initialise, air properties
    activation_and_nucleation!(work, parameters)        # 2  droplet and ice nucleation
    instant_phase_change!(work, parameters)             # 3  instant melting and freezing
    incloud_values!(work, parameters)                   # 4  in-cloud values, sizes
    autoconversion!(work, parameters)                   # 5  warm rain
    ice_snow_processes!(work, parameters)               # 6  ice, snow, vapor
    graupel_processes!(work, parameters)                # 7  graupel
    conservation_limiters!(work, parameters)            # 8  conservation
    process_tendencies!(work, parameters)               # 9  tendencies from rates
    reapply_instant_changes!(work, parameters)          # 10 stage 3 back into tendencies
    sedimentation_fall_speeds!(work, parameters)        # 11 fall speeds
    sedimentation!(work, parameters)                    # 12 sedimentation
    instant_phase_change_after_sedimentation!(work, parameters)  # 13
    effective_radius!(work, parameters)                 # 14 effective radii
    output_diagnostics!(work, parameters)               # 15 output diagnostics
    return work
end
