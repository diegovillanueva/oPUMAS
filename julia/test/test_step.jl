# 🔁🧪 M4: one full microphysics_step! against the Fortran outputs of golden_base.nc.

"""
    compare_step(work, ds) -> Vector of (name, error, column, level, identical)

🔬 Worst error of every golden variable that `output_paths()` knows: relative error,
or absolute error where the reference is zero (tolerances of handoffs/julia.plan.md),
and whether the field is bit-for-bit identical.
"""
function compare_step(work, ds)
    results = Tuple{String,Float64,Int,Int,Bool}[]
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

@testset "M4 full step against golden_base.nc" begin
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    work = allocate(read_inputs(golden_path("golden_inputs.nc")))
    microphysics_step!(work, parameters)
    NCDataset(golden_path("golden_base.nc")) do ds
        results = compare_step(work, ds)
        @test length(results) > 100
        for (name, err, column, level, _) in results
            err <= 1e-10 || println("❌ $name worst error $err at column $column level $level")
            @test err <= 1e-10
        end
        name, err, column, level, _ = results[argmax([r[2] for r in results])]
        println("🔁 worst step variable $name: $err (column $column, level $level)")
        println("🔁 bit-for-bit identical: $(count(last, results)) of $(length(results))")
    end
end
