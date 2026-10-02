# ⚖️🧪 M5 physical checks that need no Fortran: column water and energy budgets, no
#    negative mixing ratios after the step, and zero tendencies in a clear dry column.

const SPECIES = (:cloud_liquid, :cloud_ice, :rain, :snow, :graupel)

"Column integral `Σ x Δp / g` of a (column, level) field (kg m⁻² per unit of `x`)."
column_integral(x, inputs, g) = vec(sum(x .* inputs.pressure_thickness ./ g; dims = 2))

"""
    water_budget(work, parameters) -> (residual, gross)

💧 Column water budget in kg m⁻² s⁻¹: vapor plus condensate tendencies plus surface
precipitation (m s⁻¹ of liquid water, times 1000 kg m⁻³) should cancel. `gross` is the sum
of the absolute terms, the scale for rounding.
"""
function water_budget(work, parameters)
    t, inputs, g = work.tendencies, work.inputs, parameters.host.gravity
    I(x) = column_integral(x, inputs, g)
    surface = work.diagnostics.surface.total .* 1000
    residual = I(t.vapor) .+ sum(I(getproperty(t, n).q) for n in SPECIES) .+ surface
    gross = I(abs.(t.vapor)) .+ sum(I(abs.(getproperty(t, n).q)) for n in SPECIES) .+
            surface
    return residual, gross
end

"""
    energy_budget(work, parameters) -> (residual, gross)

🔥 Column energy budget in W m⁻², with ice as the reference phase: the heating `tlat`
plus `L_s` times the vapor tendency plus `L_f` times the liquid tendencies, plus `L_f`
times the liquid part of the surface precipitation that leaves the column, is zero.
"""
function energy_budget(work, parameters)
    t, inputs, h = work.tendencies, work.inputs, parameters.host
    I(x) = column_integral(x, inputs, h.gravity)
    Lv, Lf = h.latent_heat_vaporization, h.latent_heat_fusion
    surface = work.diagnostics.surface
    liquid_out = Lf .* (surface.total .- surface.ice) .* 1000
    residual = I(t.heating) .+ (Lv + Lf) .* I(t.vapor) .+
               Lf .* (I(t.cloud_liquid.q) .+ I(t.rain.q)) .+ liquid_out
    gross = I(abs.(t.heating)) .+ (Lv + Lf) .* I(abs.(t.vapor)) .+
            Lf .* (I(abs.(t.cloud_liquid.q)) .+ I(abs.(t.rain.q))) .+ abs.(liquid_out)
    return residual, gross
end

"""
    worst_negative(work) -> (excess, value, label)

➖ Most negative state after the step, `x + tendency Δt`, beyond what rounding allows:
`1e-12` of the largest `|x| + |tendency| Δt` of that variable in the column, plus
`qsmall Δt`: the conservation limiters of PUMAS only act on sink rates ≥ `qsmall`, so a
smaller sink may overdraw a tiny mixing ratio by up to `qsmall Δt` (FINDINGS.md).
`excess < 0` marks a negative value beyond that.
"""
function worst_negative(work)
    inputs, t, dt = work.inputs, work.tendencies, work.inputs.time_step
    worst = (0.0, 0.0, "")
    for n in SPECIES, m in (:q, :N)
        before = getproperty(getproperty(inputs, n), m)
        change = getproperty(getproperty(t, n), m) .* dt
        after = before .+ change
        scale = maximum(abs.(before) .+ abs.(change); dims = 2)
        for i in CartesianIndices(after)
            excess = after[i] + 1e-12 * scale[i[1]] + oPUMAS.smallest_mixing_ratio * dt
            excess < worst[1] && (worst = (excess, after[i], "$n.$m at $(Tuple(i))"))
        end
    end
    return worst
end

@testset "M5 physical checks" begin
    for label in switch_labels()
        work, parameters = stepped(label)
        water, water_gross = water_budget(work, parameters)
        energy, energy_gross = energy_budget(work, parameters)
        # 💧 qsmall-scale clipping leaves residuals of ≤ 2e-18 kg m⁻² s⁻¹ (FINDINGS.md)
        @test all(abs.(water) .<= 1e-12 .* water_gross .+ 1e-17)
        @test all(abs.(energy) .<= 1e-12 .* energy_gross .+ 1e-11)
        excess, value, where = worst_negative(work)
        @test excess >= 0
        println("⚖️ $(rpad(label, 20)) water $(maximum(abs, water)) kg m⁻² s⁻¹, energy ",
            "$(maximum(abs, energy)) W m⁻², negative beyond rounding ",
            "$excess ($value at $where)")
    end

    # 🏜️ clear dry column: no condensate, no cloud, 0.1 g/kg vapor, no host freezing
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    inputs = read_inputs(golden_path("golden_inputs.nc"))
    none() = zero(inputs.T)
    dry = merge(inputs,
        (q_vapor = fill(1e-7, size(inputs.T)),
            map(n -> n => (q = none(), N = none()), SPECIES)...,
            snow_from_host = (q = none(), N = none()),
            cloud_fraction = none(), liquid_cloud_fraction = none(),
            ice_cloud_fraction = none(),
            host_freezing = (immersion = none(), contact = none(), deposition = none())))
    work = microphysics_step!(allocate(dry), parameters)
    t = work.tendencies
    @test all(iszero, t.heating) && all(iszero, t.vapor)
    for n in SPECIES, m in (:q, :N)
        @test all(iszero, getproperty(getproperty(t, n), m))
    end
    @test all(iszero, work.diagnostics.surface.total)
end
