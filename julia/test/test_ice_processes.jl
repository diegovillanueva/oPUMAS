# ❄️🧪 M3: ice, snow, graupel, vapor and effective radius processes against the
#    captured Fortran calls.

"The rates field `key` of each NamedTuple in `rates`."
field(rates, key) = getproperty.(rates, key)

"Ambient air NamedTuple used by the vapor diffusion twins."
air_at(T, ρ, diffusivity, viscosity, schmidt, q_vapor, q_sat_liquid, q_sat_ice) =
    (; T, ρ, diffusivity, viscosity, schmidt, q_vapor, q_sat_liquid, q_sat_ice)

@testset "M3 ice, snow and graupel" begin
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    NCDataset(golden_path("golden_processes.nc")) do ds
        input(site, key) = argument(ds, site, key, :in)
        inputs(site, keys...) = (input(site, k) for k in keys)

        check_sites(ds, "ice_autoconversion") do site
            @test input(site, "a05") == parameters.settings.ice_to_snow_size
            rates = ice_to_snow_autoconversion.(Ref(parameters),
                inputs(site, "a01", "a02", "a03", "a04")...)
            [worst_error(ds, site, "a06", field(rates, :mass)),
             worst_error(ds, site, "a07", field(rates, :number))]
        end

        check_sites(ds, "snow_self_aggregation") do site
            twin = snow_self_aggregation.(Ref(parameters),
                inputs(site, "a01", "a02", "a03", "a04", "a05", "a06")...)
            [worst_error(ds, site, "a07", twin)]
        end

        check_sites(ds, "accrete_cloud_ice_snow") do site
            rates = ice_collection_by_snow.(Ref(parameters), inputs(site,
                "a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08")...)
            [worst_error(ds, site, "a09", field(rates, :mass)),
             worst_error(ds, site, "a10", field(rates, :number))]
        end

        check_sites(ds, "accrete_cloud_water_snow") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11", "a12")
            rates = droplet_collection_by_snow.(Ref(parameters), inputs(site, keys...)...)
            [worst_error(ds, site, "a13", field(rates, :mass)),
             worst_error(ds, site, "a14", field(rates, :number))]
        end

        check_sites(ds, "secondary_ice_production") do site
            rates = rime_splintering_by_snow.(inputs(site, "a01", "a02")...)
            [worst_error(ds, site, "a03", field(rates, :mass)),
             worst_error(ds, site, "a04", field(rates, :number)),
             worst_error(ds, site, "a02", field(rates, :riming))]
        end

        check_sites(ds, "accrete_rain_snow") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11", "a12")
            rates = rain_collection_by_snow.(Ref(parameters), inputs(site, keys...)...)
            [worst_error(ds, site, "a13", field(rates, :mass)),
             worst_error(ds, site, "a14", field(rates, :number))]
        end

        check_sites(ds, "heterogeneous_rain_freezing") do site
            rates = rain_immersion_freezing.(Ref(parameters),
                inputs(site, "a01", "a02", "a03", "a04")...)
            [worst_error(ds, site, "a05", field(rates, :mass)),
             worst_error(ds, site, "a06", field(rates, :number))]
        end

        check_sites(ds, "bergeron_process_snow") do site
            T, ρ, Dv, μ, Sc, qvl, qvi = inputs(site, "a01", "a02", "a03", "a04", "a05",
                "a06", "a07")
            air = air_at.(T, ρ, Dv, μ, Sc, NaN, qvl, qvi)
            twin = bergeron_to_snow.(Ref(parameters), air,
                inputs(site, "a08", "a09", "a10", "a11", "a12")...)
            [worst_error(ds, site, "a13", twin)]
        end

        check_sites(ds, "vapor_deposition_onto_snow") do site
            T, q, q_snow, N_snow, f_precip, ρ, Dv, qvl, qvi, a_snow, μ, Sc = inputs(
                site, "a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11", "a12")
            air = air_at.(T, ρ, Dv, μ, Sc, q, qvl, qvi)
            twin = snow_deposition.(Ref(parameters), air, q_snow, N_snow, f_precip,
                a_snow)
            [worst_error(ds, site, "a13", twin)]
        end

        check_sites(ds, "ice_deposition_sublimation") do site
            T, q, q_ice, N_ice, f_ice, ρ, Dv, qvl, qvi = inputs(site, "a01", "a02",
                "a03", "a04", "a05", "a06", "a07", "a08", "a09")
            air = air_at.(T, ρ, Dv, NaN, NaN, q, qvl, qvi)
            rates = ice_deposition.(Ref(parameters), air, q_ice, N_ice, f_ice)
            [worst_error(ds, site, "a10", field(rates, :bergeron)),
             worst_error(ds, site, "a11", field(rates, :deposition)),
             worst_error(ds, site, "a12", field(rates, :sublimation))]
        end

        check_sites(ds, "evaporate_sublimate_precip_graupel") do site
            @test input(site, "a31") ==
                  parameters.settings.ifs_evaporation_threshold
            T, ρ, Dv, μ, Sc, q, qvl, qvi = inputs(site, "a01", "a02", "a03", "a04",
                "a05", "a06", "a07", "a08")
            f_liquid, f_precip, a_rain, a_snow, a_graupel, b_graupel = inputs(site,
                "a09", "a10", "a11", "a12", "a13", "a14")
            q_cloud, q_ice, q_rain, q_snow, q_graupel = inputs(site, "a15", "a16",
                "a17", "a18", "a19")
            λr, N₀r, λs, N₀s, λg, N₀g = inputs(site, "a20", "a21", "a22", "a23",
                "a24", "a25")
            air = air_at.(T, ρ, Dv, μ, Sc, q, qvl, qvi)
            rain = map((q, λ, N₀, a) -> (; q, λ, N₀, a), q_rain, λr, N₀r, a_rain)
            snow = map((q, λ, N₀, a) -> (; q, λ, N₀, a), q_snow, λs, N₀s, a_snow)
            graupel = map((q, λ, N₀, a) -> (; q, λ, N₀, a, b = b_graupel),
                q_graupel, λg, N₀g, a_graupel)
            rates = precipitation_evaporation.(Ref(parameters), air, f_liquid,
                f_precip, q_cloud, q_ice, rain, snow, graupel)
            [worst_error(ds, site, "a26", field(rates, :rain)),
             worst_error(ds, site, "a27", field(rates, :snow)),
             worst_error(ds, site, "a28", field(rates, :graupel)),
             worst_error(ds, site, "a29", field(rates, :area))]
        end

        check_sites(ds, "graupel_collecting_snow") do site
            keys = ("a03", "a04", "a05", "a01", "a02", "a06", "a07", "a08", "a09")
            twin = snow_collection_by_rain.(inputs(site, keys...)...)
            [worst_error(ds, site, "a10", twin)]
        end

        check_sites(ds, "graupel_collecting_cld_water") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08")
            rates = droplet_collection_by_graupel.(inputs(site, keys...)...)
            [worst_error(ds, site, "a09", field(rates, :mass)),
             worst_error(ds, site, "a10", field(rates, :number))]
        end

        check_sites(ds, "graupel_riming_liquid_snow") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11")
            rates = rimed_snow_to_graupel.(Ref(parameters), inputs(site, keys...)...)
            [worst_error(ds, site, "a12", field(rates, :mass)),
             worst_error(ds, site, "a13", field(rates, :number)),
             worst_error(ds, site, "a01", field(rates, :riming))]
        end

        check_sites(ds, "graupel_collecting_rain") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11")
            rates = rain_collection_by_graupel.(inputs(site, keys...)...)
            [worst_error(ds, site, "a12", field(rates, :mass)),
             worst_error(ds, site, "a13", field(rates, :number))]
        end

        check_sites(ds, "graupel_rain_riming_snow") do site
            keys = ("a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08", "a09",
                "a10", "a11", "a12")
            rates = rain_snow_collisions_to_graupel.(inputs(site, keys...)...)
            [worst_error(ds, site, "a13", field(rates, :mass)),
             worst_error(ds, site, "a14", field(rates, :number)),
             worst_error(ds, site, "a01", field(rates, :rain_by_snow)),
             worst_error(ds, site, "a02", field(rates, :rain_number_by_snow)),
             worst_error(ds, site, "a03", field(rates, :snow_by_rain))]
        end

        check_sites(ds, "graupel_rime_splintering") do site
            rates = rime_splintering_by_graupel.(Ref(parameters),
                inputs(site, "a01", "a02", "a03", "a04", "a05", "a06")...)
            [worst_error(ds, site, "a07", field(rates, :droplet_mass)),
             worst_error(ds, site, "a08", field(rates, :droplet_number)),
             worst_error(ds, site, "a09", field(rates, :rain_mass)),
             worst_error(ds, site, "a10", field(rates, :rain_number)),
             worst_error(ds, site, "a05", field(rates, :droplet_riming)),
             worst_error(ds, site, "a06", field(rates, :rain_riming))]
        end

        check_sites(ds, "calc_rercld") do site
            keys = ("a08", "a01", "a02", "a03", "a04", "a05", "a06", "a07")
            twin = add_rain_cloud_effective_radius.(inputs(site, keys...)...)
            [worst_error(ds, site, "a08", twin)]
        end
    end
end
