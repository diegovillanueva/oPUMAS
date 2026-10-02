# ⏱️🧪 M8: 48 steps of the time loop against golden_loop.nc, written by
#    tests/column/driver.F90 with nsteps = 48 under the base settings.

@testset "M8 time loop against golden_loop.nc" begin
    if isfile(joinpath(GOLDEN, "golden_loop.nc"))
        parameters = setup(read_settings(golden_path("settings_base.nml")))
        work = allocate(read_inputs(golden_path("golden_inputs.nc")))
        time_loop!(work, parameters, 48)
        NCDataset(golden_path("golden_loop.nc")) do ds
            results = compare_step(work, ds)
            @test length(results) > 100
            for (name, err, column, level, _) in results
                err <= 1e-10 ||
                    println("❌ loop $name error $err at column $column level $level")
                @test err <= 1e-10
            end
            println("⏱️ loop of 48 steps   bit-for-bit $(count(last, results)) of ",
                "$(length(results))")
        end
    else
        @warn "⏱️ golden_loop.nc missing: run tests/column/run.sh"
    end
end
