# 📏🧪 M3: size distributions against the captured Fortran calls.

@testset "M3 size distributions" begin
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    NCDataset(golden_path("golden_processes.nc")) do ds
        check_sites(ds, "size_dist_param_basic", site -> begin
            h = category(ds, site, "a01", parameters)
            q = argument(ds, site, "a02", :in)
            N = argument(ds, site, "a03", :in)
            sd = size_distribution.(Ref(h), q, N)
            results = [
                worst_error(ds, site, "a03", getproperty.(sd, :number)),
                worst_error(ds, site, "a04", getproperty.(sd, :slope)),
            ]
            if haskey(ds, "$(site)__n0__out")
                push!(results, worst_error(ds, site, "n0", getproperty.(sd, :intercept)))
            end
            results
        end)

        check_sites(ds, "size_dist_param_liq", site -> begin
            h = category(ds, site, "a01", parameters)
            q = argument(ds, site, "a02", :in)
            N = argument(ds, site, "a03", :in)
            ρ = argument(ds, site, "a04", :in)
            sd = droplet_size_distribution.(Ref(h), q, N, ρ)
            [
                worst_error(ds, site, "a03", getproperty.(sd, :number)),
                worst_error(ds, site, "a05", getproperty.(sd, :shape)),
                worst_error(ds, site, "a06", getproperty.(sd, :slope)),
            ]
        end)

        check_sites(ds, "avg_diameter_vec", site -> begin
            q = argument(ds, site, "a01", :in)
            N = argument(ds, site, "a02", :in)
            ρ = argument(ds, site, "a03", :in)
            ρₚ = argument(ds, site, "a04", :in)
            [worst_error(ds, site, "a05", mean_diameter.(q, N, ρ, ρₚ))]
        end)
    end
end
