# 🔣📖 Symbols used inside oPUMAS.jl functions

🔬 Short symbols are local only; every one used in `src/` has a row here. Names that cross function boundaries are full words (see `PHILOSOPHY.md`).

| Symbol | Long name | PUMAS name | Units |
|---|---|---|---|
| `T`, `T₀` | air temperature, melting point | `t`, `tmelt` | K |
| `p` | air pressure | `p` | Pa |
| `e` | vapor pressure | `es` | Pa |
| `ε` | molar mass ratio of water vapor to dry air | `epsilo` | 1 |
| `L` | latent heat | `xxlv`, `xxls` | J kg⁻¹ |
| `q` | mixing ratio of the category in question | `qcic`, `qric`, … | kg kg⁻¹ |
| `N` | number concentration of the category in question | `ncic`, `nric`, … | kg⁻¹ |
| `ρ`, `ρ_air` | air density | `rho` | kg m⁻³ |
| `ρₚ` | bulk density of a particle | `rhow`, `rhoi`, … | kg m⁻³ |
| `λ` | slope of the gamma size distribution | `lam*` | m⁻¹ |
| `μ` | shape parameter of the droplet size distribution | `pgam` | 1 |
| `N₀` | intercept of the size distribution | `n0*` | m⁻⁴ (kg⁻¹ m⁻¹ per air mass) |
| `c`, `d` | coefficient and exponent of the mass law m = c Dᵈ | `shape_coef`, `eff_dim` | kg m⁻ᵈ, 1 |
| `a`, `b` | prefactor and exponent of the fall speed law V = a Dᵇ | `a*`, `b*` | m¹⁻ᵇ s⁻¹, 1 |
| `ν` | inverse relative variance of sub-grid cloud water | `relvar` | 1 |
| `ν` (in `sb2001_autoconversion`) | mass spectral shape of droplets | `nu` | 1 |
| `τ` | rain fraction of liquid water q_r / (q_c + q_r) | `dum` | 1 |
| `Φ` | Seifert-Beheng universal function of τ | `dum1` | 1 |
| `x_star` | drop mass separating cloud from rain | `2.6e-7` | g |
| `k` | lower index into the shape table | `dumi` | 1 |
| `U` | mass- or number-weighted fall speed | `umr`, `uns`, … | m s⁻¹ |
| `E` | collection efficiency | `eii`, `ecid`, `ecr`, `eci` | 1 |
| `D` | mean droplet diameter (μ + 1) / λ | `dc0` | m |
| `x` | autoconversion size over mean ice size λ D_cs | `d_rat` | 1 |
| `r` | slope ratio λ_r / λ_s | `d_rat` | 1 |
| `T_3`, `T_5`, `T_8` | 3, 5 and 8 K below the melting point | `tm_3`, `tm_5`, `tm_8` | K |
| `λr`, `λg` | rain and graupel slopes (local to `rain_collection_by_graupel`) | `lamr`, `lamg` | m⁻¹ |
| `ρ_s²`, `ρ_w²` | squared bulk densities of snow and water | `cons18`, `cons19` | kg² m⁻⁶ |
| `Γ` | gamma function value of the fall speed exponent | `gamma_half_*` | 1 |

## 🌍🔣 Constants and time step of the stages

| Symbol | Meaning | PUMAS name | Units |
|---|---|---|---|
| `Δt`, `inv_Δt` | time step and its inverse | `deltat`, `rdeltat` | s, s⁻¹ |
| `g` | gravity | `g` | m s⁻² |
| `Rd`, `Rv` | gas constants of dry air and water vapor | `r`, `rv` | J kg⁻¹ K⁻¹ |
| `cp` | specific heat of dry air | `cpp` | J kg⁻¹ K⁻¹ |
| `Lv`, `Lf`, `Ls` | latent heats of vaporization, fusion, sublimation | `xxlv`, `xlf`, `xxls` | J kg⁻¹ |
| `ρ_ref` | reference air density of the fall speed laws | `rhosu` | kg m⁻³ |
| `T_ice_nucleation`, `T_snow_melt`, `T_rain_freeze` | temperature thresholds | `icenuct`, `snowmelt`, `rainfrze` | K |
| `ρ_water`, `ρ_cloud_ice`, `ρ_solid_water`, `ρ_snow` | bulk densities (`stages/common.jl`) | `rhow`, `rhoi`, `rhows`, `rhosn` | kg m⁻³ |
| `smallest_mixing_ratio` | below this a category counts as absent | `qsmall` | kg kg⁻¹ |
| `one_minus_small` | number just below one, against round-off | `omsm` | 1 |
| `minimum_cloud_fraction` | floor of cloud and precipitation fractions | `mincld` | 1 |
| `new_ice_crystal_mass` | mass of a newly nucleated crystal | `mi0` | kg |
| `minimum_frozen_droplet_mass` | smallest mass of a frozen droplet | `mi0l_min` | kg |
| `number_sublimation_factor` | number loss per mass loss in sublimation | `sublim_factor` | 1 |

## 🧮🏷️ Named intermediates of the stages

🔁 Fortran scratch variables (`dum`, `dum1`, `dum2`, `dum_2D`, …) became names that say what they hold; `ratio` is kept where it is a limiter ratio. Where a name reuses a scratch array, the binding comment says so (`# dum_2D, here …`).

| Name | Stages | Holds |
|---|---|---|
| `melt_fraction`, `freeze_fraction` | s03, s13 | share of a category that melts or freezes this step |
| `latent_heating` | s03, s13 | heating rate of that phase change |
| `ΔT_all` | s03, s13 | temperature change if all of it changed phase |
| `condensation_liquid`, `condensation`, `ice_share`, `qv_after_step` | s13 | supersaturation removal |
| `sink_over_step`, `deposition_total`, `deposition_max`, `vapor_excess`, `qv_after_deposition` | s08 | vapor and mass limiters |
| `nr_sinks`, `nr_sources`, `nc_immersion_counted` | s08, s09 | rain number limiter and the counted immersion freezing |
| `rain_share` | s08 | rain share of the evaporation (reused for the ice share) |
| `qr_for_accretion`, `qc_for_accretion`, `nc_for_accretion`, `λ_power` | s06 | inputs of the rain–snow and droplet–graupel collection |
| `mass_gamma_top/bottom`, `number_gamma_top/bottom` | s11 | gamma arguments of the weighted fall speeds |
| `layer_depth`, `fall_depth`, `q_new`, `inflow_share`, `flux_divergence`, `substeps`, `steps`, `inv_steps` | s12 | sedimentation |
| `before`, `number_before` | s10, s14 | a field before a reset or limiter |
| `z_cloud_liquid`, `z_rain`, `z_cloud_ice`, `z_snow`, `z_graupel`, `z_ice_and_snow`, `z_mm6`, `condensate_total`, `ρ_dense`, `rain_rate_threshold`, `z_min_10cm` | s15 | radar reflectivity |

## 🧺📋 Stage locals and their container paths

🔗 Each stage binds its fields to readable locals (`(; inputs, state, ...) = work`, then `qc_autoconversion = rates.cloud_liquid.autoconversion_to_rain.mass  # prc`). The table is every such binding in `src/stages/`, generated from the source; the PUMAS name is the binding's comment, or the local itself when the two agree.

| Local | Container path | PUMAS name | Stages |
|---|---|---|---|
| `ni_available` | `diagnostics.cloud_ice.available_nuclei` | ncai | s02 |
| `nc_activated` | `diagnostics.cloud_liquid.activation` | ncal | s02 |
| `qc_sink_rate` | `diagnostics.cloud_liquid.first_order_sink_rate` | qcsinksum_rate1ord | s09 |
| `qc_limiter_ratio` | `diagnostics.cloud_liquid.limiter_ratio` | qcrat | s01, s08 |
| `diameter_cloud_ice` | `diagnostics.effective_diameter.cloud_ice` | deffi | s01, s14 |
| `radius_cloud_ice` | `diagnostics.effective_radius.cloud_ice` | effi | s01, s14 |
| `radius_cloud_liquid` | `diagnostics.effective_radius.cloud_liquid` | effc | s01, s14 |
| `radius_cloud_liquid_fixed` | `diagnostics.effective_radius.cloud_liquid_fixed_number` | effc_fn | s01, s14 |
| `radius_graupel` | `diagnostics.effective_radius.graupel` | reff_grau | s15 |
| `radius_rain` | `diagnostics.effective_radius.rain` | reff_rain | s15 |
| `radius_rain_and_cloud` | `diagnostics.effective_radius.rain_and_cloud` | rercld | s09 |
| `radius_snow` | `diagnostics.effective_radius.snow` | reff_snow | s15 |
| `flux_cloud_ice` | `diagnostics.flux.cloud_ice` | iflx | s12 |
| `flux_cloud_liquid` | `diagnostics.flux.cloud_liquid` | lflx | s12 |
| `flux_graupel` | `diagnostics.flux.graupel` | gflx | s12 |
| `flux_rain` | `diagnostics.flux.rain` | rflx | s12, s15 |
| `flux_snow` | `diagnostics.flux.snow` | sflx | s12 |
| `ng_per_volume` | `diagnostics.graupel.N_per_volume` | ngout | s09, s15 |
| `ng_weighted` | `diagnostics.graupel.N_weighted` | ngout2 | s15 |
| `frequency_graupel` | `diagnostics.graupel.frequency` | freqg | s15 |
| `diameter_graupel` | `diagnostics.graupel.mean_diameter` | dgout2 | s15 |
| `qg_out` | `diagnostics.graupel.q` | qgout | s09, s15 |
| `qg_weighted` | `diagnostics.graupel.q_weighted` | qgout2 | s15 |
| `diameter_graupel_radiation` | `diagnostics.graupel.radiation_diameter` | dgout | s15 |
| `ice_mass_fraction` | `diagnostics.ice_mass_fraction` | nfice | s15 |
| `precip_evaporation` | `diagnostics.precipitation.evaporation` | nevapr | s09, s10 |
| `evaporation_area` | `diagnostics.precipitation.evaporation_area` | am_evp_st | s07 |
| `precip_production` | `diagnostics.precipitation.production` | prain | s09, s10 |
| `cloudsat_dbz` | `diagnostics.radar.cloudsat_dbz` | csrfl | s15 |
| `cloudsat_dbz_weighted` | `diagnostics.radar.cloudsat_dbz_weighted` | acsrfl | s15 |
| `cloudsat_frequency` | `diagnostics.radar.cloudsat_frequency` | fcsrfl | s15 |
| `dbz_10cm` | `diagnostics.radar.reflectivity_10cm_dbz` | refl10cm | s01, s15 |
| `z_10cm` | `diagnostics.radar.reflectivity_10cm_z` | reflz10cm | s15 |
| `dbz` | `diagnostics.radar.reflectivity_dbz` | refl | s01, s15 |
| `dbz_weighted` | `diagnostics.radar.reflectivity_dbz_weighted` | arefl | s15 |
| `dbz_frequency` | `diagnostics.radar.reflectivity_frequency` | frefl | s15 |
| `z_weighted` | `diagnostics.radar.reflectivity_z_weighted` | areflz | s15 |
| `μc_radiation` | `diagnostics.radiation.droplet_shape` | pgamrad | s14 |
| `λc_radiation` | `diagnostics.radiation.droplet_slope` | lamcrad | s14 |
| `nr_per_volume` | `diagnostics.rain.N_per_volume` | nrout | s09, s14, s15 |
| `nr_weighted` | `diagnostics.rain.N_weighted` | nrout2 | s15 |
| `rain_evaporation` | `diagnostics.rain.evaporation` | prer_evap | s09 |
| `frequency_rain` | `diagnostics.rain.frequency` | freqr | s15 |
| `diameter_rain` | `diagnostics.rain.mean_diameter` | drout2 | s15 |
| `qr_out` | `diagnostics.rain.q` | qrout | s09, s15 |
| `qr_weighted` | `diagnostics.rain.q_weighted` | qrout2 | s15 |
| `ns_per_volume` | `diagnostics.snow.N_per_volume` | nsout | s09, s15 |
| `ns_weighted` | `diagnostics.snow.N_weighted` | nsout2 | s15 |
| `frequency_snow` | `diagnostics.snow.frequency` | freqs | s15 |
| `diameter_snow` | `diagnostics.snow.mean_diameter` | dsout2 | s15 |
| `qs_out` | `diagnostics.snow.q` | qsout | s09, s15 |
| `qs_weighted` | `diagnostics.snow.q_weighted` | qsout2 | s15 |
| `diameter_snow_radiation` | `diagnostics.snow.radiation_diameter` | dsout | s15 |
| `precip_surface_ice` | `diagnostics.surface.ice` | preci | s12 |
| `precip_surface` | `diagnostics.surface.total` | prect | s12 |
| `area_density_ice` | `diagnostics.surface_area_density.cloud_ice` | sadice | s14 |
| `area_density_snow` | `diagnostics.surface_area_density.snow` | sadsnow | s14 |
| `net_condensation` | `diagnostics.vapor.net_condensation` | cmeout | s09, s13 |
| `cloud_fraction` | `incloud.cloud_fraction` | cldm | s01, s04, s06 |
| `ni_cloud` | `incloud.cloud_ice.N` | niic | s04, s06 |
| `N₀i` | `incloud.cloud_ice.N₀` | n0i | s06 |
| `qi_cloud` | `incloud.cloud_ice.q` | qiic | s04, s06, s07 |
| `λi` | `incloud.cloud_ice.λ` | lami | s06, s11, s14, s15 |
| `nc_cloud` | `incloud.cloud_liquid.N` | ncic | s04, s05, s06, s07, s09 |
| `qc_cloud` | `incloud.cloud_liquid.q` | qcic | s04, s05, s06, s07, s09 |
| `λc` | `incloud.cloud_liquid.λ` | lamc | s04, s06, s09, s11, s14 |
| `μc` | `incloud.cloud_liquid.μ` | pgam | s04, s05, s06, s09, s11, s14 |
| `ng_precip` | `incloud.graupel.N` | ngic | s06 |
| `N₀g` | `incloud.graupel.N₀` | n0g | s06, s07, s15 |
| `vg_number` | `incloud.graupel.fall_speed_number` | ung | s06, s07, s11 |
| `qg_precip` | `incloud.graupel.q` | qgic | s06, s07 |
| `λg` | `incloud.graupel.λ` | lamg | s06, s07, s11, s14, s15 |
| `ice_fraction` | `incloud.ice_fraction` | icldm | s01, s02, s04, s06, s08, s09, s11, s12, s13, s14, s15 |
| `liquid_fraction` | `incloud.liquid_fraction` | lcldm | s01, s04, s06, s07, s08, s09, s11, s12, s13, s14, s15 |
| `precip_fraction` | `incloud.precipitation_fraction` | precip_frac | s01, s04, s06, s07, s08, s09, s11, s13, s14, s15 |
| `nr_precip` | `incloud.rain.N` | nric | s04, s06, s07, s09 |
| `N₀r` | `incloud.rain.N₀` | n0r | s04, s06, s07, s09, s14, s15 |
| `vr_number` | `incloud.rain.fall_speed_number` | unr | s06, s07, s11 |
| `qr_precip` | `incloud.rain.q` | qric | s04, s05, s06, s07, s09 |
| `λr` | `incloud.rain.λ` | lamr | s04, s06, s07, s08, s09, s11, s13, s14, s15 |
| `qsat_factor` | `incloud.saturation_scaling` | qsfm | s01 |
| `ns_precip` | `incloud.snow.N` | nsic | s06, s07 |
| `N₀s` | `incloud.snow.N₀` | n0s | s06, s07 |
| `vs_number` | `incloud.snow.fall_speed_number` | uns | s06, s11 |
| `qs_precip` | `incloud.snow.q` | qsic | s06, s07 |
| `λs` | `incloud.snow.λ` | lams | s06, s07, s11, s14, s15 |
| `T` | `inputs.T` | t | s01, s02, s03, s06, s07, s08, s11, s13 |
| `accretion_enhancement` | `inputs.accretion_enhancement` | accre_enhan | s06 |
| `droplet_activation` | `inputs.activated_droplets` | npccn | s02, s10 |
| `ice_nuclei` | `inputs.activated_ice_nuclei` | naai | s02 |
| `cloud_fraction_in` | `inputs.cloud_fraction` | cldn | s01 |
| `ni_in` | `inputs.cloud_ice.N` | nin | s01 |
| `qi_in` | `inputs.cloud_ice.q` | qin | s01 |
| `nc_in` | `inputs.cloud_liquid.N` | ncn | s01, s10 |
| `qc_in` | `inputs.cloud_liquid.q` | qcn | s01 |
| `ng_in` | `inputs.graupel.N` | ngr | s01, s10 |
| `qg_in` | `inputs.graupel.q` | qgr | s01, s10 |
| `contact_freezing_host` | `inputs.host_freezing.contact` | frzcnt | s06 |
| `deposition_freezing_host` | `inputs.host_freezing.deposition` | frzdep | s06 |
| `immersion_freezing_host` | `inputs.host_freezing.immersion` | frzimm | s06 |
| `ice_fraction_in` | `inputs.ice_cloud_fraction` | icecldf | s01 |
| `ice_radius_in` | `inputs.ice_effective_radius` | re_ice | s14 |
| `p_interface` | `inputs.interface_pressure` | pint | s11 |
| `liquid_fraction_in` | `inputs.liquid_cloud_fraction` | liqcldf | s01 |
| `p` | `inputs.pressure` | p | s01, s08, s13 |
| `Δp` | `inputs.pressure_thickness` | pdel | s01, s12 |
| `qv` | `inputs.q_vapor` | q | s01, s06, s08, s13 |
| `nr_in` | `inputs.rain.N` | nrn | s01, s10 |
| `qr_in` | `inputs.rain.q` | qrn | s01, s10 |
| `relative_variance` | `inputs.relative_variance` | relvar | s05, s06 |
| `qsat_factor_in` | `inputs.saturation_scaling` | qsatfac | s01 |
| `ns_in` | `inputs.snow.N` | nsn | s01, s10 |
| `qs_in` | `inputs.snow.q` | qsn | s01, s10 |
| `ns_from_host` | `inputs.snow_from_host.N` | tnd_nsnow | s06 |
| `qs_from_host` | `inputs.snow_from_host.q` | tnd_qsnow | s06 |
| `Δt` | `inputs.time_step` | deltat | s02, s06, s07, s09, s11, s12, s13, s14, s15 |
| `qi_accretion_snow_mean` | `process_rates.cloud_ice.accretion_by_snow.mass` | praitot | s09 |
| `ni_accretion_snow_mean` | `process_rates.cloud_ice.accretion_by_snow.number` | npraitot | s09 |
| `qi_autoconversion_mean` | `process_rates.cloud_ice.autoconversion_to_snow.mass` | prcitot | s09 |
| `ni_autoconversion_mean` | `process_rates.cloud_ice.autoconversion_to_snow.number` | nprcitot | s09 |
| `vi_mass` | `process_rates.cloud_ice.fall_speed.mass` | vtrmi | s11 |
| `qi_melting_mean` | `process_rates.cloud_ice.melting.mass` | melttot | s13 |
| `ni_melting_mean` | `process_rates.cloud_ice.melting.number` | nmelttot | s13 |
| `qi_net_deposition` | `process_rates.cloud_ice.net_deposition.mass` | cmeitot | s09 |
| `qi_nucleation_mean` | `process_rates.cloud_ice.nucleation_from_vapor.mass` | mnuccdtot | s09 |
| `ni_nucleation_mean` | `process_rates.cloud_ice.nucleation_from_vapor.number` | nnuccdtot | s09 |
| `qi_residual_deposition` | `process_rates.cloud_ice.residual_deposition.mass` | qirestot | s13 |
| `qi_sedimentation` | `process_rates.cloud_ice.sedimentation.mass` | qisedten | s12 |
| `ni_sedimentation` | `process_rates.cloud_ice.sedimentation.number` | nisedten | s12 |
| `qi_sedimentation_sublimation` | `process_rates.cloud_ice.sedimentation_sublimation.mass` | qisevap | s12 |
| `qi_splintering_graupel_mean` | `process_rates.cloud_ice.splintering_by_graupel.mass` | qmultgtot | s09 |
| `ni_splintering_graupel_mean` | `process_rates.cloud_ice.splintering_by_graupel.number` | nmultgtot | s09 |
| `qi_splintering_snow_mean` | `process_rates.cloud_ice.splintering_by_snow.mass` | msacwitot | s09 |
| `ni_splintering_snow_mean` | `process_rates.cloud_ice.splintering_by_snow.number` | nsacwitot | s09 |
| `qi_splintering_rain_mean` | `process_rates.cloud_ice.splintering_rain_graupel.mass` | qmultrgtot | s09 |
| `ni_splintering_rain_mean` | `process_rates.cloud_ice.splintering_rain_graupel.number` | nmultrgtot | s09 |
| `qc_accretion_mean` | `process_rates.cloud_liquid.accretion_by_rain.mass` | pratot | s09 |
| `nc_accretion_mean` | `process_rates.cloud_liquid.accretion_by_rain.number` | npratot | s09 |
| `qc_autoconversion_mean` | `process_rates.cloud_liquid.autoconversion_to_rain.mass` | prctot | s09 |
| `qc_bergeron_ice_mean` | `process_rates.cloud_liquid.bergeron_to_cloud_ice.mass` | bergtot | s09 |
| `qc_bergeron_snow_mean` | `process_rates.cloud_liquid.bergeron_to_snow.mass` | bergstot | s09 |
| `qc_contact_freezing_mean` | `process_rates.cloud_liquid.contact_freezing.mass` | mnuccttot | s09 |
| `nc_contact_freezing_mean` | `process_rates.cloud_liquid.contact_freezing.number` | nnuccttot | s09 |
| `qc_deposition_nucleation_mean` | `process_rates.cloud_liquid.deposition_nucleation.mass` | mnudeptot | s09 |
| `nc_deposition_nucleation_mean` | `process_rates.cloud_liquid.deposition_nucleation.number` | nnudeptot | s09 |
| `vc_mass` | `process_rates.cloud_liquid.fall_speed.mass` | vtrmc | s11 |
| `qc_homogeneous_freezing_mean` | `process_rates.cloud_liquid.homogeneous_freezing.mass` | homotot | s13 |
| `nc_homogeneous_freezing_mean` | `process_rates.cloud_liquid.homogeneous_freezing.number` | nhomotot | s13 |
| `qc_immersion_freezing_mean` | `process_rates.cloud_liquid.immersion_freezing.mass` | mnuccctot | s09 |
| `nc_immersion_freezing_mean` | `process_rates.cloud_liquid.immersion_freezing.number` | nnuccctot | s09 |
| `qc_residual_condensation` | `process_rates.cloud_liquid.residual_condensation.mass` | qcrestot | s13 |
| `qc_riming_graupel_mean` | `process_rates.cloud_liquid.riming_by_graupel.mass` | psacwgtot | s09 |
| `nc_riming_graupel_mean` | `process_rates.cloud_liquid.riming_by_graupel.number` | npsacwgtot | s09 |
| `qc_riming_snow_mean` | `process_rates.cloud_liquid.riming_by_snow.mass` | psacwstot | s09 |
| `nc_riming_snow_mean` | `process_rates.cloud_liquid.riming_by_snow.number` | npsacwstot | s09 |
| `qc_sedimentation` | `process_rates.cloud_liquid.sedimentation.mass` | qcsedten | s12 |
| `nc_sedimentation` | `process_rates.cloud_liquid.sedimentation.number` | ncsedten | s12 |
| `qc_sedimentation_evaporation` | `process_rates.cloud_liquid.sedimentation_evaporation.mass` | qcsevap | s12 |
| `μc_out` | `process_rates.cloud_liquid.size_distribution.shape` | pgam_out | s04 |
| `λc_out` | `process_rates.cloud_liquid.size_distribution.slope` | lamc_out | s04 |
| `qc_kk2000` | `process_rates.cloud_liquid.warm_rain_kk2000.mass` | qctend_KK2000 | s05, s06 |
| `nc_kk2000` | `process_rates.cloud_liquid.warm_rain_kk2000.number` | nctend_KK2000 | s05, s06 |
| `qc_sb2001` | `process_rates.cloud_liquid.warm_rain_sb2001.mass` | qctend_SB2001 | s05, s06 |
| `nc_sb2001` | `process_rates.cloud_liquid.warm_rain_sb2001.number` | nctend_SB2001 | s05, s06 |
| `qg_deposition_mean` | `process_rates.graupel.deposition.mass` | prdgtot | s09 |
| `vg_mass` | `process_rates.graupel.fall_speed.mass` | umg | s06, s07, s11 |
| `qg_rimed_snow_mean` | `process_rates.graupel.from_snow_riming.mass` | pgsacwtot | s09 |
| `ng_rimed_snow_mean` | `process_rates.graupel.from_snow_riming.number` | nscngtot | s09 |
| `qg_melting_mean` | `process_rates.graupel.melting.mass` | meltgtot | s03 |
| `ng_melting_mean` | `process_rates.graupel.melting.number` | nmeltgtot | s09 |
| `qg_rain_snow_mean` | `process_rates.graupel.rain_snow_collection.mass` | pgracstot | s09 |
| `ng_rain_snow_mean` | `process_rates.graupel.rain_snow_collection.number` | ngracstot | s09 |
| `qg_sedimentation` | `process_rates.graupel.sedimentation.mass` | qgsedten | s12 |
| `ng_sedimentation` | `process_rates.graupel.sedimentation.number` | ngsedten | s12 |
| `heating_rain_freezing` | `process_rates.heating.rain_freezing` | frzrdttot | s03, s13 |
| `heating_snow_melting` | `process_rates.heating.snow_melting` | meltsdttot | s03, s13 |
| `qr_accretion_graupel_mean` | `process_rates.rain.accretion_by_graupel.mass` | pracgtot | s09 |
| `nr_accretion_graupel_mean` | `process_rates.rain.accretion_by_graupel.number` | npracgtot | s09 |
| `qr_accretion_snow_mean` | `process_rates.rain.accretion_by_snow.mass` | pracstot | s09 |
| `nr_accretion_snow_mean` | `process_rates.rain.accretion_by_snow.number` | npracstot | s09 |
| `nr_autoconversion_mean` | `process_rates.rain.autoconversion_from_cloud_liquid.number` | nprctot | s09 |
| `vr_mass` | `process_rates.rain.fall_speed.mass` | umr | s06, s07, s11 |
| `qr_freezing_ice_mean` | `process_rates.rain.freezing_to_cloud_ice.mass` | mnuccritot | s09 |
| `nr_freezing_ice_mean` | `process_rates.rain.freezing_to_cloud_ice.number` | nnuccritot | s09 |
| `qr_freezing_snow_mean` | `process_rates.rain.freezing_to_snow.mass` | mnuccrtot | s09 |
| `nr_freezing_snow_mean` | `process_rates.rain.freezing_to_snow.number` | nnuccrtot | s09 |
| `qr_sedimentation` | `process_rates.rain.sedimentation.mass` | qrsedten | s12 |
| `nr_sedimentation` | `process_rates.rain.sedimentation.number` | nrsedten | s12 |
| `nr_self_collection_mean` | `process_rates.rain.self_collection.number` | nraggtot | s09 |
| `N₀r_out` | `process_rates.rain.size_distribution.intercept` | n0r_out | s04 |
| `λr_out` | `process_rates.rain.size_distribution.slope` | lamr_out | s04 |
| `qr_kk2000` | `process_rates.rain.warm_rain_kk2000.mass` | qrtend_KK2000 | s05, s06 |
| `nr_kk2000` | `process_rates.rain.warm_rain_kk2000.number` | nrtend_KK2000 | s05, s06 |
| `qr_sb2001` | `process_rates.rain.warm_rain_sb2001.mass` | qrtend_SB2001 | s05, s06 |
| `nr_sb2001` | `process_rates.rain.warm_rain_sb2001.number` | nrtend_SB2001 | s05, s06 |
| `qs_collection_rain_mean` | `process_rates.snow.collection_by_rain.mass` | psacrtot | s09 |
| `qs_deposition_mean` | `process_rates.snow.deposition.mass` | vapdepstot | s09 |
| `vs_mass` | `process_rates.snow.fall_speed.mass` | ums | s06, s07, s11 |
| `qs_melting_mean` | `process_rates.snow.melting.mass` | meltstot | s03 |
| `ns_melting_mean` | `process_rates.snow.melting.number` | nmeltstot | s09 |
| `qs_production` | `process_rates.snow.production.mass` | prodsnow | s09, s10 |
| `qs_sedimentation` | `process_rates.snow.sedimentation.mass` | qssedten | s12 |
| `ns_sedimentation` | `process_rates.snow.sedimentation.number` | nssedten | s12 |
| `qs_sublimation_mean` | `process_rates.snow.sublimation.mass` | evapsnow | s09, s10 |
| `qv_residual_condensation` | `process_rates.vapor.residual_condensation.mass` | qvres | s13 |
| `qi_accretion_snow` | `rates.cloud_ice.accretion_by_snow.mass` | prai | s06, s08, s09 |
| `ni_accretion_snow` | `rates.cloud_ice.accretion_by_snow.number` | nprai | s06, s08, s09 |
| `qi_autoconversion` | `rates.cloud_ice.autoconversion_to_snow.mass` | prci | s06, s08, s09 |
| `ni_autoconversion` | `rates.cloud_ice.autoconversion_to_snow.number` | nprci | s06, s08, s09 |
| `qi_deposition` | `rates.cloud_ice.deposition.mass` | vap_dep | s06, s08, s09 |
| `qi_nucleation` | `rates.cloud_ice.nucleation_from_vapor.mass` | mnuccd | s02, s08, s09 |
| `ni_nucleation` | `rates.cloud_ice.nucleation_from_vapor.number` | nnuccd | s02, s08, s09 |
| `qi_splintering_graupel` | `rates.cloud_ice.splintering_by_graupel.mass` | qmultg | s07, s08, s09 |
| `ni_splintering_graupel` | `rates.cloud_ice.splintering_by_graupel.number` | nmultg | s07, s08, s09 |
| `qi_splintering_snow` | `rates.cloud_ice.splintering_by_snow.mass` | msacwi | s06, s08, s09 |
| `ni_splintering_snow` | `rates.cloud_ice.splintering_by_snow.number` | nsacwi | s06, s08, s09 |
| `qi_splintering_rain` | `rates.cloud_ice.splintering_rain_graupel.mass` | qmultrg | s07, s08, s09 |
| `ni_splintering_rain` | `rates.cloud_ice.splintering_rain_graupel.number` | nmultrg | s07, s08, s09 |
| `qi_sublimation` | `rates.cloud_ice.sublimation.mass` | ice_sublim | s06, s08, s09 |
| `ni_sublimation` | `rates.cloud_ice.sublimation.number` | nsubi | s06, s08, s09 |
| `qc_accretion` | `rates.cloud_liquid.accretion_by_rain.mass` | pra | s06, s08, s09 |
| `nc_accretion` | `rates.cloud_liquid.accretion_by_rain.number` | npra | s06, s08, s09 |
| `qc_autoconversion` | `rates.cloud_liquid.autoconversion_to_rain.mass` | prc | s05, s06, s08, s09 |
| `nc_autoconversion` | `rates.cloud_liquid.autoconversion_to_rain.number` | nprc1 | s05, s08, s09 |
| `qc_bergeron_ice` | `rates.cloud_liquid.bergeron_to_cloud_ice.mass` | berg | s06, s08, s09 |
| `qc_bergeron_snow` | `rates.cloud_liquid.bergeron_to_snow.mass` | bergs | s06, s08, s09 |
| `qc_contact_freezing` | `rates.cloud_liquid.contact_freezing.mass` | mnucct | s06, s08, s09 |
| `nc_contact_freezing` | `rates.cloud_liquid.contact_freezing.number` | nnucct | s06, s08, s09 |
| `qc_deposition_nucleation` | `rates.cloud_liquid.deposition_nucleation.mass` | mnudep | s06, s08, s09 |
| `nc_deposition_nucleation` | `rates.cloud_liquid.deposition_nucleation.number` | nnudep | s06, s08, s09 |
| `nc_evaporation` | `rates.cloud_liquid.evaporation.number` | nsubc | s06, s08, s09 |
| `qc_immersion_freezing` | `rates.cloud_liquid.immersion_freezing.mass` | mnuccc | s06, s08, s09 |
| `nc_immersion_freezing` | `rates.cloud_liquid.immersion_freezing.number` | nnuccc | s06, s08, s09 |
| `qc_riming_graupel` | `rates.cloud_liquid.riming_by_graupel.mass` | psacwg | s07, s08, s09 |
| `nc_riming_graupel` | `rates.cloud_liquid.riming_by_graupel.number` | npsacwg | s07, s08, s09 |
| `qc_riming_snow` | `rates.cloud_liquid.riming_by_snow.mass` | psacws | s06, s07, s08, s09 |
| `nc_riming_snow` | `rates.cloud_liquid.riming_by_snow.number` | npsacws | s06, s08, s09 |
| `qg_deposition` | `rates.graupel.deposition.mass` | prdg | s07, s08, s09 |
| `qg_rimed_snow` | `rates.graupel.from_snow_riming.mass` | pgsacw | s07, s08, s09 |
| `ng_rimed_snow` | `rates.graupel.from_snow_riming.number` | nscng | s07, s08, s09 |
| `qg_instant_melting` | `rates.graupel.instant_melting.mass` | minstgm | s03 |
| `ng_instant_melting` | `rates.graupel.instant_melting.number` | ninstgm | s03, s09 |
| `qg_rain_snow` | `rates.graupel.rain_snow_collection.mass` | pgracs | s07, s08, s09 |
| `ng_rain_snow` | `rates.graupel.rain_snow_collection.number` | ngracs | s07, s08, s09 |
| `qr_accretion_graupel` | `rates.rain.accretion_by_graupel.mass` | pracg | s07, s08, s09 |
| `nr_accretion_graupel` | `rates.rain.accretion_by_graupel.number` | npracg | s07, s08, s09 |
| `qr_accretion_snow` | `rates.rain.accretion_by_snow.mass` | pracs | s06, s07, s08, s09 |
| `nr_accretion_snow` | `rates.rain.accretion_by_snow.number` | npracs | s06, s07, s08, s09 |
| `nr_autoconversion` | `rates.rain.autoconversion_from_cloud_liquid.number` | nprc | s05, s06, s08, s09 |
| `qr_evaporation` | `rates.rain.evaporation.mass` | pre | s07, s08, s09 |
| `nr_evaporation` | `rates.rain.evaporation.number` | nsubr | s08, s09 |
| `qr_freezing_ice` | `rates.rain.freezing_to_cloud_ice.mass` | mnuccri | s08, s09 |
| `nr_freezing_ice` | `rates.rain.freezing_to_cloud_ice.number` | nnuccri | s08, s09 |
| `qr_freezing_snow` | `rates.rain.freezing_to_snow.mass` | mnuccr | s06, s08, s09 |
| `nr_freezing_snow` | `rates.rain.freezing_to_snow.number` | nnuccr | s06, s08, s09 |
| `qr_instant_freezing` | `rates.rain.instant_freezing.mass` | minstrf | s03 |
| `nr_instant_freezing` | `rates.rain.instant_freezing.number` | ninstrf | s03 |
| `nr_self_collection` | `rates.rain.self_collection.number` | nragg | s06, s08, s09 |
| `ns_aggregation` | `rates.snow.aggregation.number` | nsagg | s06, s08, s09 |
| `qs_collection_rain` | `rates.snow.collection_by_rain.mass` | psacr | s07, s08, s09 |
| `qs_deposition` | `rates.snow.deposition.mass` | vap_deps | s06, s08, s09 |
| `qs_instant_melting` | `rates.snow.instant_melting.mass` | minstsm | s03 |
| `ns_instant_melting` | `rates.snow.instant_melting.number` | ninstsm | s03, s09 |
| `qs_sublimation` | `rates.snow.sublimation.mass` | prds | s07, s08, s09 |
| `ns_sublimation` | `rates.snow.sublimation.number` | nsubs | s08, s09 |
| `graupel_share` | `scratch.evaporation_share.graupel` | dum3A | s08 |
| `rain_share` | `scratch.evaporation_share.rain` | dum1A | s08 |
| `snow_share` | `scratch.evaporation_share.snow` | dum2A | s08 |
| `ni_intercept` | `scratch.ice_intercept` | dumni0A2D | s14, s15 |
| `z_interface` | `scratch.interface_height` | zint | s11, s12 |
| `inv_Δp` | `scratch.inverse_pressure_thickness` | pdel_inv | s01, s12 |
| `λ_power` | `scratch.number_before` | dum_2D, here λ to the fall speed exponent | s06 |
| `qv_after_deposition` | `scratch.number_before` | dum_2D, here the vapor after deposition | s08 |
| `before` | `scratch.number_before` | dum_2D, here the value before the reset | s10 |
| `qv_after_step` | `scratch.number_before` | dum_2D, here the vapor after the step | s13 |
| `number_before` | `scratch.number_before` | dum_2D | s14 |
| `condensate_total` | `scratch.number_before` | dum_2D, here the total condensate | s15 |
| `fall_ni` | `scratch.pressure_fall_speed.cloud_ice.N` | fni | s11, s12 |
| `fall_qi` | `scratch.pressure_fall_speed.cloud_ice.q` | fi | s11, s12 |
| `fall_nc` | `scratch.pressure_fall_speed.cloud_liquid.N` | fnc | s11, s12 |
| `fall_qc` | `scratch.pressure_fall_speed.cloud_liquid.q` | fc | s11, s12 |
| `fall_ng` | `scratch.pressure_fall_speed.graupel.N` | fng | s11, s12 |
| `fall_qg` | `scratch.pressure_fall_speed.graupel.q` | fg | s11, s12 |
| `fall_nr` | `scratch.pressure_fall_speed.rain.N` | fnr | s11, s12 |
| `fall_qr` | `scratch.pressure_fall_speed.rain.q` | fr | s11, s12 |
| `fall_ns` | `scratch.pressure_fall_speed.snow.N` | fns | s11, s12 |
| `fall_qs` | `scratch.pressure_fall_speed.snow.q` | fs | s11, s12 |
| `qsat_adjusted` | `scratch.q_sat_adjusted` | qvnA | s08, s13 |
| `qsat_estimate` | `scratch.q_sat_estimate` | qvnAI | s08 |
| `heating_sed_ice` | `scratch.sedimentation_heating.cloud_ice` | tlat_i | s12 |
| `heating_sed_liquid` | `scratch.sedimentation_heating.cloud_liquid` | tlat_l | s12 |
| `qv_sed_ice` | `scratch.sedimentation_vapor.cloud_ice` | qvlat_i | s12 |
| `qv_sed_liquid` | `scratch.sedimentation_vapor.cloud_liquid` | qvlat_l | s12 |
| `ni_sed` | `scratch.sedimenting.cloud_ice.N` | dumni | s11, s12, s13, s14 |
| `qi_sed` | `scratch.sedimenting.cloud_ice.q` | dumi | s11, s12, s13, s14, s15 |
| `nc_sed` | `scratch.sedimenting.cloud_liquid.N` | dumnc | s11, s12, s13, s14 |
| `qc_sed` | `scratch.sedimenting.cloud_liquid.q` | dumc | s11, s12, s13, s14 |
| `ng_sed` | `scratch.sedimenting.graupel.N` | dumng | s11, s12, s13, s14 |
| `qg_sed` | `scratch.sedimenting.graupel.q` | dumg | s11, s12, s13, s14 |
| `nr_sed` | `scratch.sedimenting.rain.N` | dumnr | s11, s12, s13, s14 |
| `qr_sed` | `scratch.sedimenting.rain.q` | dumr | s11, s12, s13, s14 |
| `ns_sed` | `scratch.sedimenting.snow.N` | dumns | s11, s12, s13, s14 |
| `qs_sed` | `scratch.sedimenting.snow.q` | dums | s11, s12, s13, s14 |
| `ns_intercept` | `scratch.snow_intercept` | dumns0A2D | s14, s15 |
| `surface_frozen_cloud_ice` | `scratch.surface_ice.cloud_ice` | preci_i | s12 |
| `surface_frozen_graupel` | `scratch.surface_ice.graupel` | preci_g | s12 |
| `surface_frozen_snow` | `scratch.surface_ice.snow` | preci_s | s12 |
| `surface_cloud_ice` | `scratch.surface_precipitation.cloud_ice` | prect_i | s12 |
| `surface_cloud_liquid` | `scratch.surface_precipitation.cloud_liquid` | prect_l | s12 |
| `surface_graupel` | `scratch.surface_precipitation.graupel` | prect_g | s12 |
| `surface_rain` | `scratch.surface_precipitation.rain` | prect_r | s12 |
| `surface_snow` | `scratch.surface_precipitation.snow` | prect_s | s12 |
| `T_estimate` | `scratch.temperature_estimate` | ttmpA | s08, s13 |
| `qv_estimate` | `scratch.vapor_estimate` | qtmpAI | s08 |
| `ρ` | `state.air.density` | rho | s01, s02, s04, s05, s06, s07, s09, s11, s12, s13, s14, s15 |
| `ρ_correction` | `state.air.density_correction` | rhof | s01, s06, s11 |
| `D_v` | `state.air.diffusivity` | dv | s01 |
| `a_cloud_ice` | `state.air.fall_speed_prefactor.cloud_ice` | ain | s01, s11 |
| `a_cloud_liquid` | `state.air.fall_speed_prefactor.cloud_liquid` | acn | s01, s11 |
| `a_graupel` | `state.air.fall_speed_prefactor.graupel` | agn | s01, s06, s07, s11 |
| `a_rain` | `state.air.fall_speed_prefactor.rain` | arn | s01, s06, s07, s11 |
| `a_small_ice` | `state.air.fall_speed_prefactor.small_cloud_ice` | ajn | s01, s11 |
| `a_snow` | `state.air.fall_speed_prefactor.snow` | asn | s01, s06, s07, s11 |
| `qsat_ice` | `state.air.q_sat_ice` | qvi | s01, s08 |
| `qsat_liquid` | `state.air.q_sat_liquid` | qvl | s01 |
| `rh` | `state.air.relative_humidity` | relhum | s01, s02 |
| `Sc` | `state.air.schmidt` | sc | s01 |
| `e_sat_ice` | `state.air.vapor_pressure_ice` | esi | s01, s02 |
| `e_sat_liquid` | `state.air.vapor_pressure_liquid` | esl | s01, s02 |
| `η` | `state.air.viscosity` | mu | s01, s06 |
| `ni` | `state.cloud_ice.N` | ni | s01, s04, s06, s08, s11, s13, s14 |
| `qi` | `state.cloud_ice.q` | qi | s01, s04, s06, s08, s11, s13, s14, s15 |
| `nc` | `state.cloud_liquid.N` | nc | s01, s02, s04, s08, s09, s10, s11, s13, s14, s15 |
| `qc` | `state.cloud_liquid.q` | qc | s01, s02, s04, s08, s09, s11, s13, s14, s15 |
| `ng` | `state.graupel.N` | ng | s01, s03, s06, s09, s10, s11, s13, s14 |
| `qg` | `state.graupel.q` | qg | s01, s03, s04, s06, s08, s09, s10, s11, s13, s14 |
| `nr` | `state.rain.N` | nr | s01, s03, s04, s08, s09, s10, s11, s13, s14 |
| `qr` | `state.rain.q` | qr | s01, s03, s04, s08, s09, s10, s11, s13, s14 |
| `ns` | `state.snow.N` | ns | s01, s03, s06, s08, s09, s10, s11, s13, s14 |
| `qs` | `state.snow.q` | qs | s01, s03, s04, s06, s08, s09, s10, s11, s13, s14 |
| `ni_tendency` | `tendencies.cloud_ice.N` | nitend | s09, s11, s12, s13, s14 |
| `qi_tendency` | `tendencies.cloud_ice.q` | qitend | s09, s11, s12, s13, s14 |
| `nc_tendency` | `tendencies.cloud_liquid.N` | nctend | s09, s10, s11, s12, s13, s14, s15 |
| `qc_tendency` | `tendencies.cloud_liquid.q` | qctend | s09, s11, s12, s13, s14 |
| `ng_tendency` | `tendencies.graupel.N` | ngtend | s09, s10, s11, s12, s13, s14 |
| `qg_tendency` | `tendencies.graupel.q` | qgtend | s09, s10, s11, s12, s13, s14 |
| `heating` | `tendencies.heating` | tlat | s03, s09, s12, s13 |
| `nr_tendency` | `tendencies.rain.N` | nrtend | s09, s10, s11, s12, s13, s14 |
| `qr_tendency` | `tendencies.rain.q` | qrtend | s09, s10, s11, s12, s13, s14 |
| `ns_tendency` | `tendencies.snow.N` | nstend | s09, s10, s11, s12, s13, s14 |
| `qs_tendency` | `tendencies.snow.q` | qstend | s09, s10, s11, s12, s13, s14 |
| `qv_tendency` | `tendencies.vapor` | qvlat | s09, s12, s13 |
