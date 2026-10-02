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
