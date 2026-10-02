# 🔬🧷 Helpers to replay the process calls captured inside micro_pumas_tend.

"Names of the captured call sites of `routine` (`routine_1`, `routine_2`, ...)."
function call_sites(ds, routine)
    sites = String[]
    for n in 1:99
        site = "$(routine)_$(n)"
        any(startswith("$(site)__"), keys(ds)) || break
        push!(sites, site)
    end
    return sites
end

"One argument of a call site, before (`:in`) or after (`:out`) the call."
argument(ds, site, key, phase) = captured(ds, "$(site)__$(key)__$(phase)")

"The hydrometeor category whose properties were passed as argument `key`."
function category(ds, site, key, parameters)
    expr = ds["$(site)__$(key)_rho__in"].attrib["expr"]
    name = Dict("mg_liq_props%rho" => :cloud_liquid, "mg_ice_props%rho" => :cloud_ice,
        "mg_rain_props%rho" => :rain, "mg_snow_props%rho" => :snow,
        "mg_graupel_props%rho" => :graupel, "mg_hail_props%rho" => :hail)[expr]
    return parameters.hydrometeors[name]
end

"""
Worst relative error of `twin` against the Fortran output `key` of a site; exact zeros
in the Fortran must be exact zeros in the twin.
"""
function worst_error(ds, site, key, twin)
    fortran = argument(ds, site, key, :out)
    errors = relative_error.(twin, fortran)
    worst, at = findmax(errors)
    return (; worst, site, key, at, twin = twin[at], fortran = fortran[at])
end

"Worst relative error per checked routine, printed at the end of the tests."
const REPORT = Pair{String,Float64}[]

"Collects the worst error per routine and reports the overall worst."
function check_sites(ds, routine, replay; tolerance = 1e-13)
    sites = call_sites(ds, routine)
    isempty(sites) && return @test_skip "no captured calls of $(routine)"
    results = [r for site in sites for r in replay(site)]
    worst = results[argmax([r.worst for r in results])]
    push!(REPORT, routine => worst.worst)
    worst.worst <= tolerance || @info "worst site of $(routine)" worst
    @test worst.worst <= tolerance
    return worst
end

# 🧩📝 Lets a test write `check_sites(ds, routine) do site ... end`.
check_sites(replay::Function, ds, routine; kwargs...) =
    check_sites(ds, routine, replay; kwargs...)
