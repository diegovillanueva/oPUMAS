# 🌧️🧪 M3: warm rain processes against the captured Fortran calls.

@testset "M3 warm rain" begin
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    NCDataset(golden_path("golden_processes.nc")) do ds
        input(site, key) = argument(ds, site, key, :in)

        check_sites(ds,
            "kk2000_liq_autoconversion",
            site -> begin
                rates = kk2000_autoconversion.(Ref(parameters),
                    input(site, "a02"), input(site, "a03"), input(site, "a04"),
                    input(site, "a05"))
                [
                    worst_error(ds, site, "a06", getproperty.(rates, :mass)),
                    worst_error(ds, site, "a07", getproperty.(rates, :rain_number)),
                    worst_error(ds, site, "a08", getproperty.(rates, :droplet_number))
                ]
            end)

        check_sites(ds,
            "accrete_cloud_water_rain",
            site -> begin
                rates = accretion.(Ref(parameters), input(site, "a02"), input(site, "a03"),
                    input(site, "a04"), input(site, "a05"), input(site, "a06"))
                [
                    worst_error(ds, site, "a07", getproperty.(rates, :mass)),
                    worst_error(ds, site, "a08", getproperty.(rates, :droplet_number))
                ]
            end)

        check_sites(ds,
            "self_collection_rain",
            site -> begin
                twin = rain_self_collection.(input(site, "a01"), input(site, "a02"),
                    input(site, "a03"))
                [worst_error(ds, site, "a04", twin)]
            end)
    end

    # 🌧️🔀 Seifert-Beheng routines run only with warm_rain = sb2001.
    NCDataset(golden_path("golden_processes_sb2001.nc")) do ds
        input(site, key) = argument(ds, site, key, :in)

        # 🧮🎯 Real exponents (x^4.0) take libm pow in Fortran but power by squaring in
        #    Julia, so allow a little more than the default.
        check_sites(ds, "sb2001v2_liq_autoconversion"; tolerance = 1e-12) do site
            rates = sb2001_autoconversion.(input(site, "a01"), input(site, "a02"),
                input(site, "a03"), input(site, "a04"), input(site, "a05"))
            [
                worst_error(ds, site, "a07", getproperty.(rates, :mass)),
                worst_error(ds, site, "a08", getproperty.(rates, :rain_number)),
                worst_error(ds, site, "a09", getproperty.(rates, :droplet_number))
            ]
        end

        check_sites(ds,
            "sb2001v2_accre_cld_water_rain",
            site -> begin
                rates = sb2001_accretion.(input(site, "a01"), input(site, "a02"),
                    input(site, "a03"), input(site, "a04"))
                [
                    worst_error(ds, site, "a06", getproperty.(rates, :mass)),
                    worst_error(ds, site, "a07", getproperty.(rates, :droplet_number))
                ]
            end)
    end
end
