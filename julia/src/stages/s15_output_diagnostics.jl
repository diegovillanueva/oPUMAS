# 📡🌧️ Stage 15, output diagnostics (micro_pumas_v1.F90 lines 4307–4547): mean
#    diameters and effective radii of precipitation, radar reflectivities and the ice
#    mass fraction.

"""
    output_diagnostics!(work, parameters)

📡 Stage 15. Fills the `*out2`, `d*out*`, `freq*`, `reff_*`, radar and `nfice`
diagnostics from the final state.
"""
function output_diagnostics!(work, parameters)
    settings = parameters.settings
    ρ_dense = dense_precipitation(parameters).density
    (; state, incloud, tendencies, diagnostics, scratch) = work
    qc = state.cloud_liquid.q
    nc = state.cloud_liquid.N
    qi = state.cloud_ice.q
    ρ = state.air.density  # rho
    liquid_fraction = incloud.liquid_fraction  # lcldm
    ice_fraction = incloud.ice_fraction  # icldm
    precip_fraction = incloud.precipitation_fraction  # precip_frac
    nc_tendency = tendencies.cloud_liquid.N  # nctend
    flux_rain = diagnostics.flux.rain  # rflx
    qr_out = diagnostics.rain.q  # qrout
    nr_per_volume = diagnostics.rain.N_per_volume  # nrout
    qs_out = diagnostics.snow.q  # qsout
    ns_per_volume = diagnostics.snow.N_per_volume  # nsout
    qg_out = diagnostics.graupel.q  # qgout
    ng_per_volume = diagnostics.graupel.N_per_volume  # ngout
    qr_weighted = diagnostics.rain.q_weighted  # qrout2
    nr_weighted = diagnostics.rain.N_weighted  # nrout2
    qs_weighted = diagnostics.snow.q_weighted  # qsout2
    ns_weighted = diagnostics.snow.N_weighted  # nsout2
    qg_weighted = diagnostics.graupel.q_weighted  # qgout2
    ng_weighted = diagnostics.graupel.N_weighted  # ngout2
    diameter_rain = diagnostics.rain.mean_diameter  # drout2
    diameter_snow = diagnostics.snow.mean_diameter  # dsout2
    diameter_graupel = diagnostics.graupel.mean_diameter  # dgout2
    diameter_snow_radiation = diagnostics.snow.radiation_diameter  # dsout
    diameter_graupel_radiation = diagnostics.graupel.radiation_diameter  # dgout
    frequency_rain = diagnostics.rain.frequency  # freqr
    frequency_snow = diagnostics.snow.frequency  # freqs
    frequency_graupel = diagnostics.graupel.frequency  # freqg
    radius_rain = diagnostics.effective_radius.rain  # reff_rain
    radius_snow = diagnostics.effective_radius.snow  # reff_snow
    radius_graupel = diagnostics.effective_radius.graupel  # reff_grau
    dbz = diagnostics.radar.reflectivity_dbz  # refl
    dbz_weighted = diagnostics.radar.reflectivity_dbz_weighted  # arefl
    z_weighted = diagnostics.radar.reflectivity_z_weighted  # areflz
    dbz_frequency = diagnostics.radar.reflectivity_frequency  # frefl
    cloudsat_dbz = diagnostics.radar.cloudsat_dbz  # csrfl
    cloudsat_dbz_weighted = diagnostics.radar.cloudsat_dbz_weighted  # acsrfl
    cloudsat_frequency = diagnostics.radar.cloudsat_frequency  # fcsrfl
    dbz_10cm = diagnostics.radar.reflectivity_10cm_dbz  # refl10cm
    z_10cm = diagnostics.radar.reflectivity_10cm_z  # reflz10cm
    ice_mass_fraction = diagnostics.ice_mass_fraction  # nfice
    λr = incloud.rain.λ  # lamr
    N₀r = incloud.rain.N₀  # n0r
    λi = incloud.cloud_ice.λ  # lami
    λs = incloud.snow.λ  # lams
    λg = incloud.graupel.λ  # lamg
    N₀g = incloud.graupel.N₀  # n0g
    ni_intercept = scratch.ice_intercept  # dumni0A2D
    ns_intercept = scratch.snow_intercept  # dumns0A2D
    condensate_total = scratch.number_before  # dum_2D, here the total condensate
    qi_sed = scratch.sedimenting.cloud_ice.q  # dumi
    Δt = work.inputs.time_step
    columns, levels = grid_size(work)

    # 🌧️ rain
    for k in 1:levels, i in 1:columns
        diameter_rain[i, k] = mean_diameter(qr_out[i, k], nr_per_volume[i, k], ρ[i, k], ρ_water)
    end
    for k in 1:levels, i in 1:columns
        if qr_out[i, k] > 1.0e-7 && nr_per_volume[i, k] > 0.0
            qr_weighted[i, k] = qr_out[i, k] * precip_fraction[i, k]
            nr_weighted[i, k] = nr_per_volume[i, k] * precip_fraction[i, k]
            frequency_rain[i, k] = precip_fraction[i, k]
            radius_rain[i, k] = 1.5 * diameter_rain[i, k] * 1.0e6
        else
            qr_weighted[i, k] = 0.0
            nr_weighted[i, k] = 0.0
            diameter_rain[i, k] = 0.0
            frequency_rain[i, k] = 0.0
            radius_rain[i, k] = 0.0
        end
        if ns_per_volume[i, k] == 0.0
            ns_per_volume[i, k] = 1.0e-34
        end
    end

    # ❄️ snow
    for k in 1:levels, i in 1:columns
        diameter_snow[i, k] = mean_diameter(qs_out[i, k], ns_per_volume[i, k], ρ[i, k], ρ_snow)
    end
    for k in 1:levels, i in 1:columns
        if qs_out[i, k] > 1.0e-7 && ns_per_volume[i, k] > 0.0
            qs_weighted[i, k] = qs_out[i, k] * precip_fraction[i, k]
            ns_weighted[i, k] = ns_per_volume[i, k] * precip_fraction[i, k]
            frequency_snow[i, k] = precip_fraction[i, k]
            diameter_snow_radiation[i, k] = 3.0 * ρ_snow / ρ_solid_water * diameter_snow[i, k]
            radius_snow[i, k] = 1.5 * diameter_snow[i, k] * 1.0e6
        else
            diameter_snow_radiation[i, k] = 0.0
            qs_weighted[i, k] = 0.0
            ns_weighted[i, k] = 0.0
            diameter_snow[i, k] = 0.0
            frequency_snow[i, k] = 0.0
            radius_snow[i, k] = 0.0
        end
        if ng_per_volume[i, k] == 0.0
            ng_per_volume[i, k] = 1.0e-34
        end
    end

    # 🟤 graupel or hail
    dense_on = settings.do_hail || settings.do_graupel
    for k in 1:levels, i in 1:columns
        diameter_graupel[i, k] = dense_on ?
                       mean_diameter(qg_out[i, k], ng_per_volume[i, k], ρ[i, k], ρ_dense) : 0.0
    end
    for k in 1:levels, i in 1:columns
        if qg_out[i, k] > 1.0e-7 && ng_per_volume[i, k] > 0.0
            qg_weighted[i, k] = qg_out[i, k] * precip_fraction[i, k]
            ng_weighted[i, k] = ng_per_volume[i, k] * precip_fraction[i, k]
            frequency_graupel[i, k] = precip_fraction[i, k]
            diameter_graupel_radiation[i, k] = 3.0 * ρ_dense / ρ_solid_water * diameter_graupel[i, k]
            radius_graupel[i, k] = 1.5 * diameter_graupel[i, k] * 1.0e6
        else
            diameter_graupel_radiation[i, k] = 0.0
            qg_weighted[i, k] = 0.0
            ng_weighted[i, k] = 0.0
            diameter_graupel[i, k] = 0.0
            frequency_graupel[i, k] = 0.0
            radius_graupel[i, k] = 0.0
        end
    end

    # 📡 analytic radar reflectivity (Shupe), in-precipitation, mm⁶ m⁻³
    rain_rate_threshold = 0.0001 / 3600.0
    for k in 1:levels, i in 1:columns
        if qc[i, k] >= smallest_mixing_ratio && (nc[i, k] + nc_tendency[i, k] * Δt) > 10.0
            x = qc[i, k] / liquid_fraction[i, k] * ρ[i, k] * 1000.0
            z_cloud_liquid = (x * x) /
                  (0.109 * (nc[i, k] + nc_tendency[i, k] * Δt) / liquid_fraction[i, k] * ρ[i, k] /
                   1.0e6) * liquid_fraction[i, k] / precip_fraction[i, k]
        else
            z_cloud_liquid = 0.0
        end
        if qi[i, k] >= smallest_mixing_ratio
            z_ice_and_snow = (qi[i, k] * ρ[i, k] / ice_fraction[i, k] * 1000.0 / 0.1)^(1.0 / 0.63) *
                   ice_fraction[i, k] / precip_fraction[i, k]
        else
            z_ice_and_snow = 0.0
        end
        if qs_out[i, k] >= smallest_mixing_ratio
            z_ice_and_snow = z_ice_and_snow + (qs_out[i, k] * ρ[i, k] * 1000.0 / 0.1)^(1.0 / 0.63)
        end
        dbz[i, k] = z_cloud_liquid + z_ice_and_snow
        # 🌧️ Z = a R^b, R in mm h⁻¹ (Comstock et al. 2004)
        if flux_rain[i, k+1] >= rain_rate_threshold
            z_rain = 32.0 * (flux_rain[i, k+1] * 3600.0)^1.4
        else
            z_rain = 0.0
        end
        dbz[i, k] = dbz[i, k] + z_rain
        z_weighted[i, k] = dbz[i, k] * precip_fraction[i, k]
        if dbz[i, k] > radar.minimum_reflectivity
            dbz[i, k] = 10.0 * log10(dbz[i, k])
        else
            dbz[i, k] = -9999.0
        end
        if dbz[i, k] > radar.minimum_dbz
            dbz_weighted[i, k] = dbz[i, k] * precip_fraction[i, k]
            dbz_frequency[i, k] = precip_fraction[i, k]
        else
            dbz_weighted[i, k] = 0.0
            z_weighted[i, k] = 0.0
            dbz_frequency[i, k] = 0.0
        end
        cloudsat_dbz[i, k] = min(radar.cloudsat_max_dbz, dbz[i, k])
        if cloudsat_dbz[i, k] > radar.cloudsat_min_dbz
            cloudsat_dbz_weighted[i, k] = dbz[i, k] * precip_fraction[i, k]
            cloudsat_frequency[i, k] = precip_fraction[i, k]
        else
            cloudsat_dbz_weighted[i, k] = 0.0
            cloudsat_frequency[i, k] = 0.0
        end
    end

    # 📡 10 cm reflectivity (Morrison; ice dielectric factor from Smith 1984)
    z_min_10cm = radar.minimum_reflectivity_10cm
    for k in 1:levels, i in 1:columns
        z_rain = z_min_10cm
        z_cloud_ice = z_min_10cm
        z_snow = z_min_10cm
        z_graupel = z_min_10cm
        if λr[i, k] > 0.0
            x = λr[i, k]
            z_rain = ρ[i, k] * N₀r[i, k] * 720.0 / (x * x * x) / (x * x * x) / x
            z_rain = max(z_rain, z_min_10cm)
        end
        if λi[i, k] > 0.0
            x = λi[i, k]
            x2 = x * x
            z_cloud_ice = ρ[i, k] * (0.176 / 0.93) * 720.0 * ni_intercept[i, k] *
                   ((ρ_cloud_ice / 900.0) * (ρ_cloud_ice / 900.0)) / ((x * x2) * (x2 * x2))
            z_cloud_ice = max(z_cloud_ice, z_min_10cm)
        end
        if λs[i, k] > 0.0
            z_snow = ρ[i, k] * (0.176 / 0.93) * 720.0 * ns_intercept[i, k] *
                   ((ρ_snow / 900.0) * (ρ_snow / 900.0)) / λs[i, k]^7.0
            z_snow = max(z_snow, z_min_10cm)
        end
        # 🟤 `do_hail .or. do_graupel .and. lamg > 0` binds as hail or (graupel and λ > 0)
        if settings.do_hail || (settings.do_graupel && λg[i, k] > 0.0)
            z_graupel = ρ[i, k] * (0.176 / 0.93) * 720.0 * N₀g[i, k] *
                   ((ρ_dense / 900.0) * (ρ_dense / 900.0)) / λg[i, k]^7.0
            z_graupel = max(z_graupel, z_min_10cm)
        end
        z_10cm[i, k] = (z_rain + z_cloud_ice + z_snow + z_graupel) * precip_fraction[i, k]
        z_mm6 = z_10cm[i, k] * 1.0e18
        dbz_10cm[i, k] = 10.0 * log10(z_mm6)

        # 🧊 ice mass fraction
        condensate_total[i, k] = qs_out[i, k] + qr_out[i, k] + qc[i, k] + qi[i, k]
        qi_sed[i, k] = qs_out[i, k] + qi[i, k]
        if qi_sed[i, k] > smallest_mixing_ratio && condensate_total[i, k] > smallest_mixing_ratio
            ice_mass_fraction[i, k] = min(qi_sed[i, k] / condensate_total[i, k], 1.0)
        else
            ice_mass_fraction[i, k] = 0.0
        end
    end
    return nothing
end
