# 🧱🔢 M2: constants, gamma, saturation and setup against the Fortran module state.

using oPUMAS
using oPUMAS: bulk_density, fall_speed_law, ventilation, collection_efficiency
using oPUMAS: bigg_freezing, radar, fall_speed_factor

const LIMITER_OFF = 0x7FF1111111111111

"A Fortran limit marked off (`no_limiter`) is `nothing` in the twin."
fortran_limit(x::Float64) = reinterpret(UInt64, x) == LIMITER_OFF ? nothing : x

function hydrometeor_pairs(name, h::HydrometeorProperties)
    bounds = h.slope_bounds === nothing ? (nothing, nothing) : h.slope_bounds
    return [
        "utils_mg_$(name)_props_rho" => h.density,
        "utils_mg_$(name)_props_eff_dim" => h.mass_exponent,
        "utils_mg_$(name)_props_shape_coef" => h.mass_coefficient,
        "utils_mg_$(name)_props_lambda_bounds" => bounds,
        "utils_mg_$(name)_props_min_mean_mass" => h.minimum_mean_mass
    ]
end

"Fortran module variable => twin value, for every constant `micro_pumas_init` stores."
function module_state_pairs(p::Parameters)
    s, host, hyd, γ = p.settings, p.host, p.hydrometeors, p.gamma
    law = fall_speed_law
    Lᵥ, Lₛ = host.latent_heat_vaporization, p.latent_heat_sublimation
    state = Any[
        "utils_pi" => Float64(π),
        "utils_omsm" => oPUMAS.one_minus_small,
        "utils_qsmall" => oPUMAS.smallest_mixing_ratio,
        "utils_mincld" => oPUMAS.minimum_cloud_fraction,
        "utils_icsmall" => oPUMAS.small_incloud_mixing_ratio,
        "utils_rhow" => bulk_density.liquid_water,
        "utils_rhoi" => bulk_density.cloud_ice,
        "utils_rhows" => bulk_density.solid_water,
        "utils_rhosn" => bulk_density.snow,
        "utils_rhog" => bulk_density.graupel,
        "utils_rhoh" => bulk_density.hail,
        "utils_ac" => law.cloud_liquid.a, "utils_bc" => law.cloud_liquid.b,
        "utils_ai" => law.cloud_ice.a, "utils_bi" => law.cloud_ice.b,
        "utils_aj" => law.small_cloud_ice.a, "utils_bj" => law.small_cloud_ice.b,
        "utils_ar" => law.rain.a, "utils_br" => law.rain.b,
        "utils_as" => law.snow.a, "utils_bs" => law.snow.b,
        "utils_ag" => law.graupel.a, "utils_bg" => law.graupel.b,
        "utils_ah" => law.hail.a, "utils_bh" => law.hail.b,
        "utils_mi0" => oPUMAS.new_ice_crystal_mass,
        "utils_mg0" => oPUMAS.new_graupel_mass,
        "utils_mmult" => oPUMAS.splinter_mass,
        "utils_dsph" => oPUMAS.sphere_mass_exponent,
        "utils_lam_bnd_rain" => oPUMAS.rain_slope_bounds,
        "utils_lam_bnd_snow" => oPUMAS.snow_slope_bounds,
        "utils_min_mean_mass_liq" => oPUMAS.minimum_mean_mass,
        "utils_min_mean_mass_ice" => oPUMAS.minimum_mean_mass,
        "utils_f1s" => ventilation.snow.f1, "utils_f2s" => ventilation.snow.f2,
        "utils_f1r" => ventilation.rain.f1, "utils_f2r" => ventilation.rain.f2,
        "utils_eii" => collection_efficiency.ice_snow,
        "utils_ecid" => collection_efficiency.ice_droplet,
        "utils_ecr" => collection_efficiency.droplet_rain,
        "utils_bimm" => bigg_freezing.b, "utils_aimm" => bigg_freezing.a,
        "utils_droplet_mass_25um" => oPUMAS.autoconversion_drop_mass_25um,
        "utils_droplet_mass_40um" => oPUMAS.autoconversion_drop_mass_40um,
        "utils_rv" => host.gas_constant_water_vapor,
        "utils_cpp" => host.specific_heat_dry_air,
        "utils_tmelt" => host.melting_point,
        "utils_ra" => host.gas_constant_dry_air,
        "utils_xxlv" => Lᵥ,
        "utils_xlf" => host.latent_heat_fusion,
        "utils_xxls" => Lₛ,
        "utils_gamma_bs_plus3" => γ.snow_plus3,
        "utils_gamma_half_br_plus5" => γ.half_rain_plus5,
        "utils_gamma_half_bs_plus5" => γ.half_snow_plus5,
        "utils_gamma_2bs_plus2" => γ.twice_snow_plus2,
        "v1_nccons" => s.constant_droplet_number,
        "v1_nicons" => s.constant_ice_number,
        "v1_ngcons" => s.constant_graupel_number,
        "v1_nrcons" => false, "v1_nscons" => false,
        "v1_ncnst" => s.droplet_number,
        "v1_ninst" => s.ice_number,
        "v1_ngnst" => s.graupel_number,
        "v1_nrnst" => 0.0, "v1_nsnst" => 0.0,
        "v1_evap_sed_off" => s.no_sedimentation_evaporation,
        "v1_icenuc_rh_off" => s.ice_nucleation_without_rh,
        "v1_icenuc_use_meyers" => s.meyers_ice_nucleation,
        "v1_evap_scl_ifs" => s.ifs_evaporation_scaling,
        "v1_evap_rhthrsh_ifs" => s.ifs_evaporation_threshold,
        "v1_rainfreeze_ifs" => s.ifs_rain_freezing,
        "v1_ifs_sed" => s.ifs_snow_fall_speed,
        "v1_precip_fall_corr" => s.precipitation_fall_correction,
        "v1_csmin" => radar.cloudsat_min_dbz,
        "v1_csmax" => radar.cloudsat_max_dbz,
        "v1_mindbz" => radar.minimum_dbz,
        "v1_minrefl" => radar.minimum_reflectivity,
        "v1_minrefl10" => radar.minimum_reflectivity_10cm,
        "v1_dcs" => s.ice_to_snow_size,
        "v1_mi0l_min" => oPUMAS.minimum_frozen_droplet_mass,
        "v1_sublim_factor" => oPUMAS.number_sublimation_factor,
        "v1_g" => host.gravity,
        "v1_r" => host.gas_constant_dry_air,
        "v1_rv" => host.gas_constant_water_vapor,
        "v1_cpp" => host.specific_heat_dry_air,
        "v1_tmelt" => host.melting_point,
        "v1_xxlv" => Lᵥ,
        "v1_xlf" => host.latent_heat_fusion,
        "v1_xxls" => Lₛ,
        "v1_rhmini" => s.minimum_rh_for_ice_cloud,
        "v1_microp_uniform" => s.uniform_subcolumns,
        "v1_do_cldice" => s.prognostic_cloud_ice,
        "v1_use_hetfrz_classnuc" => s.host_heterogeneous_freezing,
        "v1_do_hail" => s.do_hail,
        "v1_do_graupel" => s.do_graupel,
        "v1_rhosu" => p.reference_air_density,
        "v1_icenuct" => p.ice_nucleation_temperature,
        "v1_snowmelt" => p.snow_melting_temperature,
        "v1_rainfrze" => p.rain_freezing_temperature,
        "v1_gamma_br_plus1" => γ.rain_plus1, "v1_gamma_br_plus4" => γ.rain_plus4,
        "v1_gamma_bs_plus1" => γ.snow_plus1, "v1_gamma_bs_plus4" => γ.snow_plus4,
        "v1_gamma_bi_plus1" => γ.ice_plus1, "v1_gamma_bi_plus4" => γ.ice_plus4,
        "v1_gamma_bj_plus1" => γ.small_ice_plus1, "v1_gamma_bj_plus4" => γ.small_ice_plus4,
        "v1_gamma_bg_plus1" => γ.graupel_plus1, "v1_gamma_bg_plus4" => γ.graupel_plus4,
        "v1_xxlv_squared" => Lᵥ^2,
        "v1_xxls_squared" => Lₛ^2,
        "v1_micro_mg_berg_eff_factor" => s.bergeron_efficiency,
        "v1_micro_mg_accre_enhan_fact" => s.accretion_enhancement,
        "v1_micro_mg_autocon_fact" => s.autoconversion_prefactor,
        "v1_micro_mg_autocon_nd_exp" => s.autoconversion_number_exponent,
        "v1_micro_mg_autocon_lwp_exp" => s.autoconversion_water_exponent,
        "v1_micro_mg_homog_size" => s.homogeneous_freezing_size,
        "v1_micro_mg_vtrmi_factor" => s.ice_fall_speed_factor,
        "v1_micro_mg_vtrms_factor" => s.snow_fall_speed_factor,
        "v1_micro_mg_effi_factor" => s.ice_effective_radius_factor,
        "v1_micro_mg_iaccr_factor" => s.ice_accretion_factor,
        "v1_micro_mg_max_nicons" => s.maximum_ice_number,
        "v1_remove_supersat" => s.remove_supersaturation,
        "v1_vfactor" => fall_speed_factor.precipitation,
        "v1_vfac_drop" => fall_speed_factor.cloud_liquid,
        "v1_vfac_ice" => fall_speed_factor.cloud_ice,
        "v1_do_implicit_fall" => s.implicit_sedimentation,
        "v1_accre_sees_auto" => s.accretion_sees_autoconversion
    ]
    for (name, h) in pairs(hyd)
        fortran = name == :cloud_liquid ? "liq" : name == :cloud_ice ? "ice" : string(name)
        append!(state, hydrometeor_pairs(fortran, h))
    end
    return state
end

"Bitwise equality of a captured Fortran value and a twin value."
same_bits(fortran::Real, twin::Bool) = (fortran != 0) == twin
same_bits(fortran::Float64, twin) = fortran_limit(fortran) === twin
same_bits(fortran::AbstractVector, twin) = all(same_bits.(fortran, collect(twin)))

@testset "M2 foundations" begin
    parameters = setup(read_settings(golden_path("settings_base.nml")))
    NCDataset(golden_path("golden_processes.nc")) do ds
        @testset "module state equals the Fortran bit for bit" begin
            state = module_state_pairs(parameters)
            names = ["module_state__$(name)__in" for (name, _) in state]
            fortran_names = filter(n -> startswith(n, "module_state__"), keys(ds))
            @test isempty(setdiff(fortran_names, names))   # nothing left untested
            for ((name, twin), key) in zip(state, names)
                ok = same_bits(captured(ds, key), twin)
                ok || @info "differs" name fortran=captured(ds, key) twin
                @test ok
            end
        end

        @testset "saturation matches the Fortran" begin
            sites = [(k, f)
                     for (k, f) in (("qsat_water", saturation_over_water),
                ("qsat_ice", saturation_over_ice)) for _ in 1:1]
            worst = 0.0
            for (prefix, saturation) in sites, n in 1:9
                site = "$(prefix)_$(n)"
                haskey(ds, "$(site)__a01__in") || continue
                T = captured(ds, "$(site)__a01__in")
                p = captured(ds, "$(site)__a02__in")
                e = captured(ds, "$(site)__a03__out")
                q = captured(ds, "$(site)__a04__out")
                twin = saturation.(Ref(parameters.saturation), T, p)
                err_e = maximum(relative_error.(getproperty.(twin, :vapor_pressure), e))
                err_q = maximum(relative_error.(getproperty.(twin, :specific_humidity), q))
                worst = max(worst, err_e, err_q)
            end
            @info "saturation: worst relative error" worst
            @test worst <= 1e-14
        end

        @testset "gamma is the C library's tgamma" begin
            @test gamma_function(4.0) === 6.0
            @test gamma_function(3.41) === parameters.gamma.snow_plus3
        end
    end
end
