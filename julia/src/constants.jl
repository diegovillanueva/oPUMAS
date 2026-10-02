# 🧱🔢 Fixed constants of PUMAS (module parameters of micro_pumas_utils.F90 and
#      micro_pumas_v1.F90). Settings that a namelist can change live in parameters.jl.

# 🔢 Numerical thresholds
"Smallest mixing ratio (kg kg⁻¹) that the microphysics treats as non-zero (`qsmall`)."
const smallest_mixing_ratio = 1.0e-18

"Threshold for some in-cloud mixing ratios (kg kg⁻¹) (`icsmall`)."
const small_incloud_mixing_ratio = 1.0e-8

"Number just below one, used against round-off (`omsm`)."
const one_minus_small = 1.0 - 1.0e-5

"Smallest cloud fraction used to compute in-cloud values (`mincld`)."
const minimum_cloud_fraction = 0.0001

# 🧊 Bulk densities (kg m⁻³)
"Bulk density of each hydrometeor (`rhow`, `rhoi`, `rhows`, `rhosn`, `rhog`, `rhoh`)."
const bulk_density = (
    liquid_water = 1000.0,  # rhow
    cloud_ice = 500.0,      # rhoi
    solid_water = 917.0,    # rhows
    snow = 250.0,           # rhosn
    graupel = 500.0,        # rhog
    hail = 900.0           # rhoh
)

"""
Exponent `d` of the mass-diameter law m = c Dᵈ; all particles are spheres (`dsph`).
"""
const sphere_mass_exponent = 3.0

# 🪂 Fall speed laws V = a Dᵇ (V in m s⁻¹, D in m)
"""
Prefactor `a` and exponent `b` of the fall speed law V = a Dᵇ for each hydrometeor
(`ac/bc`, `ai/bi`, `aj/bj`, `ar/br`, `as/bs`, `ag/bg`, `ah/bh`).
`small_cloud_ice` is ice below 10 μm radius, treated as a droplet of ice density.
"""
const fall_speed_law = (
    cloud_liquid = (a = 3.0e7, b = 2.0),
    cloud_ice = (a = 700.0, b = 1.0),
    # 🧊🔢 aj = ac*((rhoi/rhows)**(bc/3))*rhows/rhow, folded by the Fortran compiler.
    small_cloud_ice = (
        a = 3.0e7 * ((500.0 / 917.0)^(2.0 / 3.0)) * 917.0 / 1000.0,
        b = 2.0
    ),
    rain = (a = 841.99667, b = 0.8),
    snow = (a = 11.72, b = 0.41),
    graupel = (a = 19.3, b = 0.37),
    hail = (a = 114.5, b = 0.5)
)

# ⚖️ Particle masses (kg)
"Mass of a new ice crystal from aerosol freezing, a 1 μm ice sphere (`mi0`)."
const new_ice_crystal_mass = 4.0 / 3.0 * π * bulk_density.cloud_ice * (1.0e-6)^3

"Mass of a new graupel particle, as in Morrison et al. (2005) (`mg0`)."
const new_graupel_mass = 1.6e-10

"Mass of a 5 μm ice splinter from rime splintering (`mmult`)."
const splinter_mass = 4.0 / 3.0 * π * bulk_density.cloud_ice * (5.0e-6)^3

"Smallest mass of a newly frozen droplet, a 4 μm water sphere (`mi0l_min`)."
const minimum_frozen_droplet_mass = 4.0 / 3.0 * π * bulk_density.liquid_water * (4.0e-6)^3

"""
Mass of a raindrop made by autoconversion, 25 μm radius (`droplet_mass_25um`):
4/3 π ρ_w (25 μm)³. Written as the value the Fortran compiler folds; computing it at run
time gives one unit in the last place less (see FINDINGS.md).
"""
const autoconversion_drop_mass_25um = 6.544984694978737e-11

"Mass of a raindrop made by autoconversion, 40 μm radius (`droplet_mass_40um`)."
const autoconversion_drop_mass_40um = 4.0 / 3.0 * π * bulk_density.liquid_water *
                                      (40.0e-6)^3

# 📏 Size limits
"Bounds of the rain slope parameter λ (m⁻¹): mean diameters 500 μm to 20 μm."
const rain_slope_bounds = (1.0 / 500.0e-6, 1.0 / 20.0e-6)

"Bounds of the snow, graupel and hail slope parameter λ (m⁻¹): 2000 μm to 10 μm."
const snow_slope_bounds = (1.0 / 2000.0e-6, 1.0 / 10.0e-6)

"Smallest mean particle mass of cloud liquid and cloud ice (kg) (`min_mean_mass_*`)."
const minimum_mean_mass = 1.0e-20

# 🌬️ Ventilation coefficients f₁ + f₂ √Re Sc^(1/3)
"Ventilation coefficients of snow (`f1s`, `f2s`) and rain (`f1r`, `f2r`)."
const ventilation = (snow = (f1 = 0.86, f2 = 0.28), rain = (f1 = 0.78, f2 = 0.308))

# 🤝 Collection efficiencies
"Collection efficiencies: ice-snow aggregation, ice-droplet, droplet/rain-snow/rain."
const collection_efficiency = (
    ice_snow = 0.5,         # eii
    ice_droplet = 0.7,      # ecid
    droplet_rain = 1.0     # ecr
)

# ❄️ Immersion freezing of Bigg (1953)
"Bigg (1953) immersion freezing coefficients `bimm` (m⁻³ s⁻¹) and `aimm` (K⁻¹)."
const bigg_freezing = (b = 100.0, a = 0.66)

# 📡 Radar reflectivity
"Clamps of the radar diagnostics: CloudSat range, minimum dBZ and minimum Z."
const radar = (
    cloudsat_min_dbz = -30.0,     # csmin
    cloudsat_max_dbz = 26.0,      # csmax
    minimum_dbz = -99.0,          # mindbz
    minimum_reflectivity = 1.26e-10,   # minrefl
    minimum_reflectivity_10cm = 1.0e-26  # minrefl10
)

"Number sublimation factor (`sublim_factor`), zero: sublimation removes no particles."
const number_sublimation_factor = 0.0

"Fall speed scaling of rain, snow and graupel (`vfactor`), droplets and ice."
const fall_speed_factor = (precipitation = 1.0, cloud_liquid = 1.0, cloud_ice = 1.0)
