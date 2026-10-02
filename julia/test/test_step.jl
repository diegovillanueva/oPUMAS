# 🔁🧪 M4/M5: one full microphysics_step! against the Fortran outputs of every switch row
#    of tests/switches.csv that has a golden file (tau and wsub_default have none).

"""
    compare_step(work, ds) -> Vector of (name, error, column, level, identical)

🔬 Worst error of every golden variable that `output_paths()` knows: relative error,
or absolute error where the reference is zero (tolerances of handoffs/julia.plan.md),
and whether the field is bit-for-bit identical.
"""
function compare_step(work, ds)
    results = Tuple{String, Float64, Int, Int, Bool}[]
    for (name, (container, path...)) in output_paths()
        haskey(ds, String(name)) || continue
        v = ds[String(name)].var                     # raw: keep NaN, the fill value
        reference = v[ntuple(_ -> :, ndims(v))...]
        twin = field_at(work[container], path)
        worst, where = 0.0, (0, 0)
        for index in CartesianIndices(reference)
            a, b = twin[index], reference[index]
            err = b == 0 ? (abs(a) < 1e-20 ? 0.0 : Inf) : abs(a - b) / abs(b)
            isnan(a) && (err = isnan(b) ? 0.0 : Inf)
            if err > worst
                worst, where = err, (Tuple(index)..., 0)[1:2]
            end
        end
        push!(results, (String(name), worst, where..., isequal(twin, reference)))
    end
    return results
end

"""
    switch_labels() -> Vector{String}

🎛️ Labels of tests/switches.csv whose golden file and settings exist, in table order.
"""
function switch_labels()
    tests = joinpath(@__DIR__, "..", "..", "tests")
    rows = [readlines(joinpath(tests, "switches.csv"))[2:end];
            readlines(joinpath(tests, "column", "switches_sppe.csv"))[2:end]]
    labels = [first(split(row, ',')) for row in rows if !isempty(strip(row))]
    return filter(
        l -> isfile(joinpath(GOLDEN, "golden_$l.nc")) &&
            isfile(joinpath(GOLDEN, "settings_$l.nml")),
        labels)
end

"🪜 Workspace after one step of the golden inputs under the settings of switch row `label`."
function stepped(label)
    parameters = setup(read_settings(golden_path("settings_$label.nml")))
    work = allocate(read_inputs(golden_path("golden_inputs.nc")))
    return microphysics_step!(work, parameters), parameters
end

@testset "M5 full step against golden_<label>.nc" begin
    labels = switch_labels()
    @test "base" in labels && length(labels) >= 16
    for label in labels
        work, _ = stepped(label)
        NCDataset(golden_path("golden_$label.nc")) do ds
            results = compare_step(work, ds)
            @test length(results) > 100
            for (name, err, column, level, _) in results
                err <= 1e-10 ||
                    println("❌ $label $name error $err at column $column level $level")
                @test err <= 1e-10
            end
            name, err, _, _, _ = results[argmax([r[2] for r in results])]
            println("🔁 $(rpad(label, 20)) bit-for-bit $(count(last, results)) of ",
                "$(length(results)), worst $name: $err")
        end
    end
end
