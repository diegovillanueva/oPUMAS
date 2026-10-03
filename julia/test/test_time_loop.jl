# ⏱️🧪 M8: the time loop against golden_loop.nc, written by tests/column/driver.F90 under
#    the base settings; the step count (432 steps of 600 s = 3 days) is its attribute nsteps.

@testset "M8 time loop against golden_loop.nc" begin
    if isfile(joinpath(GOLDEN, "golden_loop.nc"))
        parameters = setup(read_settings(golden_path("settings_base.nml")))
        work = allocate(read_inputs(golden_path("golden_inputs.nc")))
        NCDataset(golden_path("golden_loop.nc")) do ds
            nsteps = Int(get(ds.attrib, "nsteps", 48))
            time_loop!(work, parameters, nsteps)
            results = compare_step(work, ds)
            @test length(results) > 100
            for (name, err, column, level, _) in results
                err <= 1e-10 ||
                    println("❌ loop $name error $err at column $column level $level")
                @test err <= 1e-10
            end
            println("⏱️ loop of $nsteps steps   bit-for-bit $(count(last, results)) of ",
                "$(length(results))")
        end
    else
        @warn "⏱️ golden_loop.nc missing: run tests/column/run.sh"
    end
end
