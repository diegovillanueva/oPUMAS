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
    work = (; inputs, state, incloud, rates, process_rates, tendencies, diagnostics,
        scratch)
    return (; work..., fortran = fortran_aliases(work))
end

"Set every array of a nested NamedTuple to zero."
zero_fields!(tree::NamedTuple) = foreach(zero_fields!, values(tree))
zero_fields!(array::AbstractArray) = fill!(array, 0.0)

"""
    fortran_aliases(work)

🏷️🧱 The matrices of `work` under their PUMAS names (`qc`, `prc`, `lamr`, ...), so that
the M4 stages read like micro_pumas_v1.F90. Every alias is the same array as its
container field. The M6 readability pass replaces these names by the container paths.
"""
function fortran_aliases(work)
    (; inputs, state, incloud, scratch) = work
    air, prefactor = state.air, state.air.fall_speed_prefactor
    aliases = Pair{Symbol,Any}[
        :t => inputs.T, :q => inputs.q_vapor,
        :qcn => inputs.cloud_liquid.q, :ncn => inputs.cloud_liquid.N,
        :qin => inputs.cloud_ice.q, :nin => inputs.cloud_ice.N,
        :qrn => inputs.rain.q, :nrn => inputs.rain.N,
        :qsn => inputs.snow.q, :nsn => inputs.snow.N,
        :qgr => inputs.graupel.q, :ngr => inputs.graupel.N,
        :relvar => inputs.relative_variance, :accre_enhan => inputs.accretion_enhancement,
        :p => inputs.pressure, :pdel => inputs.pressure_thickness,
        :pint => inputs.interface_pressure, :cldn => inputs.cloud_fraction,
        :liqcldf => inputs.liquid_cloud_fraction, :icecldf => inputs.ice_cloud_fraction,
        :qsatfac => inputs.saturation_scaling, :naai => inputs.activated_ice_nuclei,
        :npccn => inputs.activated_droplets, :tnd_qsnow => inputs.snow_from_host.q,
        :tnd_nsnow => inputs.snow_from_host.N, :re_ice => inputs.ice_effective_radius,
        :frzimm => inputs.host_freezing.immersion, :frzcnt => inputs.host_freezing.contact,
        :frzdep => inputs.host_freezing.deposition,
        :qc => state.cloud_liquid.q, :nc => state.cloud_liquid.N,
        :qi => state.cloud_ice.q, :ni => state.cloud_ice.N,
        :qr => state.rain.q, :nr => state.rain.N, :qs => state.snow.q, :ns => state.snow.N,
        :qg => state.graupel.q, :ng => state.graupel.N,
        :rho => air.density, :dv => air.diffusivity, :mu => air.viscosity,
        :sc => air.schmidt, :rhof => air.density_correction,
        :arn => prefactor.rain, :asn => prefactor.snow, :agn => prefactor.graupel,
        :acn => prefactor.cloud_liquid, :ain => prefactor.cloud_ice,
        :ajn => prefactor.small_cloud_ice,
        :esl => air.vapor_pressure_liquid, :esi => air.vapor_pressure_ice,
        :qvl => air.q_sat_liquid, :qvi => air.q_sat_ice, :relhum => air.relative_humidity,
        :qcic => incloud.cloud_liquid.q, :ncic => incloud.cloud_liquid.N,
        :pgam => incloud.cloud_liquid.μ, :lamc => incloud.cloud_liquid.λ,
        :qiic => incloud.cloud_ice.q, :niic => incloud.cloud_ice.N,
        :lami => incloud.cloud_ice.λ, :n0i => incloud.cloud_ice.N₀,
        :qric => incloud.rain.q, :nric => incloud.rain.N, :lamr => incloud.rain.λ,
        :n0r => incloud.rain.N₀, :unr => incloud.rain.fall_speed_number,
        :qsic => incloud.snow.q, :nsic => incloud.snow.N, :lams => incloud.snow.λ,
        :n0s => incloud.snow.N₀, :uns => incloud.snow.fall_speed_number,
        :qgic => incloud.graupel.q, :ngic => incloud.graupel.N, :lamg => incloud.graupel.λ,
        :n0g => incloud.graupel.N₀, :ung => incloud.graupel.fall_speed_number,
        :precip_frac => incloud.precipitation_fraction, :cldm => incloud.cloud_fraction,
        :lcldm => incloud.liquid_fraction, :icldm => incloud.ice_fraction,
        :qsfm => incloud.saturation_scaling,
        :dum_2D => scratch.number_before, :dumni0A2D => scratch.ice_intercept,
        :dumns0A2D => scratch.snow_intercept, :zint => scratch.interface_height,
        :pdel_inv => scratch.inverse_pressure_thickness,
        :tlat_i => scratch.sedimentation_heating.cloud_ice,
        :tlat_l => scratch.sedimentation_heating.cloud_liquid,
        :qvlat_i => scratch.sedimentation_vapor.cloud_ice,
        :qvlat_l => scratch.sedimentation_vapor.cloud_liquid,
        :prect_l => scratch.surface_precipitation.cloud_liquid,
        :prect_i => scratch.surface_precipitation.cloud_ice,
        :prect_r => scratch.surface_precipitation.rain,
        :prect_s => scratch.surface_precipitation.snow,
        :prect_g => scratch.surface_precipitation.graupel,
        :preci_i => scratch.surface_ice.cloud_ice, :preci_s => scratch.surface_ice.snow,
        :preci_g => scratch.surface_ice.graupel,
        :ttmpA => scratch.temperature_estimate, :qtmpAI => scratch.vapor_estimate,
        :qvnAI => scratch.q_sat_estimate, :qvnA => scratch.q_sat_adjusted,
        :dum1A => scratch.evaporation_share.rain, :dum2A => scratch.evaporation_share.snow,
        :dum3A => scratch.evaporation_share.graupel,
    ]
    for (letter, species) in (("c", :cloud_liquid), ("i", :cloud_ice), ("r", :rain),
        ("s", :snow), ("g", :graupel))
        push!(aliases, Symbol("dum", letter) => scratch.sedimenting[species].q,
            Symbol("dumn", letter) => scratch.sedimenting[species].N,
            Symbol("f", letter) => scratch.pressure_fall_speed[species].q,
            Symbol("fn", letter) => scratch.pressure_fall_speed[species].N)
    end
    for (table, container) in ((LOCAL_RATE_PATHS, work.rates),
        (PROCESS_RATE_PATHS, work.process_rates), (TENDENCY_PATHS, work.tendencies),
        (DIAGNOSTIC_PATHS, work.diagnostics), (FLUX_PATHS, work.diagnostics),
        (SURFACE_PATHS, work.diagnostics))
        for (name, path) in pairs(table)
            push!(aliases, name => field_at(container, path))
        end
    end
    names = first.(aliases)
    allunique(names) || error("duplicate PUMAS alias names")
    return NamedTuple{Tuple(names)}(Tuple(last.(aliases)))
end

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
