# 🧺🧱 The containers of one time step (julia/DESIGN.md): nested NamedTuples of
#    (column, level) matrices, allocated once and overwritten by every step.

"""
    nested_fields(paths, make)

Nested NamedTuple with one leaf `make()` per path (a tuple of Symbols), in table order.
"""
function nested_fields(paths, make)
    heads = unique(first.(paths))
    values = map(heads) do head
        rests = [path[2:end] for path in paths if first(path) == head]
        any(isempty, rests) ? make() : nested_fields(rests, make)
    end
    return NamedTuple{Tuple(heads)}(Tuple(values))
end

"Mass `q` and number `N` matrices of each prognostic hydrometeor."
function species_fields(make)
    pair() = (q = make(), N = make())
    return (cloud_liquid = pair(), cloud_ice = pair(), rain = pair(), snow = pair(),
        graupel = pair())
end

"""
    ColumnInputs(; time_step, fields...)

Inputs from the host, never modified (`intent(in)` of micro_pumas_tend). Keywords are
the PUMAS argument names: `t, q, qcn, ncn, qin, nin, qrn, nrn, qsn, nsn, qgr, ngr,
relvar, accre_enhan, p, pdel, pint, cldn, liqcldf, icecldf, qsatfac, naai, npccn,
tnd_qsnow, tnd_nsnow, re_ice, frzimm, frzcnt, frzdep`.
"""
function ColumnInputs(; time_step, f...)
    return (
        time_step,                                      # deltatin
        T = f[:t],
        q_vapor = f[:q],
        cloud_liquid = (q = f[:qcn], N = f[:ncn]),
        cloud_ice = (q = f[:qin], N = f[:nin]),
        rain = (q = f[:qrn], N = f[:nrn]),
        snow = (q = f[:qsn], N = f[:nsn]),
        graupel = (q = f[:qgr], N = f[:ngr]),
        relative_variance = f[:relvar],
        accretion_enhancement = f[:accre_enhan],
        pressure = f[:p],
        pressure_thickness = f[:pdel],
        interface_pressure = f[:pint],                  # (column, level + 1)
        cloud_fraction = f[:cldn],
        liquid_cloud_fraction = f[:liqcldf],
        ice_cloud_fraction = f[:icecldf],
        saturation_scaling = f[:qsatfac],
        activated_ice_nuclei = f[:naai],                # kg⁻¹ s⁻¹
        activated_droplets = f[:npccn],                 # kg⁻¹ s⁻¹
        snow_from_host = (q = f[:tnd_qsnow], N = f[:tnd_nsnow]),
        ice_effective_radius = f[:re_ice],              # m
        host_freezing = (immersion = f[:frzimm], contact = f[:frzcnt],
            deposition = f[:frzdep]),                   # cm⁻³ s⁻¹
    )
end

"""
    allocate(inputs) -> Workspace

The seven containers of DESIGN.md plus `scratch` (the Fortran work arrays that several
stages share), sized like `inputs.T`.
"""
function allocate(inputs)
    columns, levels = size(inputs.T)
    level() = zeros(columns, levels)
    interface() = zeros(columns, levels + 1)
    surface() = zeros(columns)
    air = (
        density = level(), diffusivity = level(), viscosity = level(),
        schmidt = level(), density_correction = level(),       # rho dv mu sc rhof
        fall_speed_prefactor = (rain = level(), snow = level(), graupel = level(),
            cloud_liquid = level(), cloud_ice = level(), small_cloud_ice = level()),
        vapor_pressure_liquid = level(), vapor_pressure_ice = level(),  # esl esi
        q_sat_liquid = level(), q_sat_ice = level(),                    # qvl qvi
        relative_humidity = level(),                                    # relhum
    )
    state = (; species_fields(level)..., air)
    distribution() = (q = level(), N = level(), λ = level(), N₀ = level(),
        fall_speed_number = level())
    incloud = (
        cloud_liquid = (q = level(), N = level(), μ = level(), λ = level()),
        cloud_ice = distribution(), rain = distribution(), snow = distribution(),
        graupel = distribution(),
        precipitation_fraction = level(), cloud_fraction = level(),
        liquid_fraction = level(), ice_fraction = level(), saturation_scaling = level(),
    )
    rates = nested_fields(collect(values(LOCAL_RATE_PATHS)), level)
    process_rates = nested_fields(collect(values(PROCESS_RATE_PATHS)), level)
    tendencies = nested_fields(collect(values(TENDENCY_PATHS)), level)
    diagnostics = (
        nested_fields(collect(values(DIAGNOSTIC_PATHS)), level)...,
        flux = nested_fields([p[2:end] for p in values(FLUX_PATHS)], interface),
        surface = nested_fields([p[2:end] for p in values(SURFACE_PATHS)], surface),
    )
    scratch = (
        sedimenting = species_fields(level),           # dumc, dumnc, ..., dumng
        pressure_fall_speed = species_fields(level),   # fc, fnc, ..., fng
        number_before = level(),                       # dum_2D
        ice_intercept = level(), snow_intercept = level(),  # dumni0A2D dumns0A2D
        interface_height = interface(),                # zint
        inverse_pressure_thickness = level(),          # pdel_inv
        sedimentation_heating = (cloud_ice = level(), cloud_liquid = level()),
        sedimentation_vapor = (cloud_ice = level(), cloud_liquid = level()),
        surface_precipitation = (cloud_liquid = surface(), cloud_ice = surface(),
            rain = surface(), snow = surface(), graupel = surface()),
        surface_ice = (cloud_ice = surface(), snow = surface(), graupel = surface()),
        temperature_estimate = level(),                # ttmpA
        vapor_estimate = level(),                      # qtmpAI
        q_sat_estimate = level(),                      # qvnAI
        q_sat_adjusted = level(),                      # qvnA
        evaporation_share = (rain = level(), snow = level(), graupel = level()),  # dum1A..3A
    )
    return (; inputs, state, incloud, rates, process_rates, tendencies, diagnostics,
        scratch)
end

"Set every array of a nested NamedTuple to zero."
zero_fields!(tree::NamedTuple) = foreach(zero_fields!, values(tree))
zero_fields!(array::AbstractArray) = fill!(array, 0.0)

"""
    read_inputs(path) -> ColumnInputs

📥🌍 Inputs of one call of micro_pumas_tend from a file written by tests/column/run.sh
(`golden_inputs.nc`): arrays become (column, level), the time step is the `deltat`
global attribute.
"""
function read_inputs(path)
    NCDataset(path) do ds
        names = (:t, :q, :qcn, :ncn, :qin, :nin, :qrn, :nrn, :qsn, :nsn, :qgr, :ngr,
            :relvar, :accre_enhan, :p, :pdel, :pint, :cldn, :liqcldf, :icecldf, :qsatfac,
            :naai, :npccn, :tnd_qsnow, :tnd_nsnow, :re_ice, :frzimm, :frzcnt, :frzdep)
        # 🧾 raw values: NaN is both the fill value and a real input value
        fields = NamedTuple{names}(map(n -> Array{Float64}(ds[String(n)].var[:, :]), names))
        return ColumnInputs(; time_step = Float64(ds.attrib["deltat"]), fields...)
    end
end
