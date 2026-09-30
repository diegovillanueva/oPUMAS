# ❄️🔬 PUMAS immersion freezing scheme

📍📂 Where the immersion freezing of cloud droplets is implemented in oPUMAS (branch `dev`, as of 2026-09-30).

## 🧮❄️ Bigg (1953) parameterization (default, `use_hetfrz_classnuc = .false.`)

- 🧩 Subroutine `immersion_freezing` in `micro_pumas_utils.F90` (~line 1581).
- 🎛️ Parameters (`micro_pumas_utils.F90` ~line 235):
  ```
  bimm = 100.0
  aimm = 0.66
  ```
- 🌡️ Active only when `qcic >= qsmall` and `t < 269.15 K` (-4 C).
- 📐 Tendencies (gamma droplet size distribution with shape `pgam`, slope `lamc`):
  ```
  nnuccc = pi/6 * ncic * (pgam+1)(pgam+2)(pgam+3) * bimm*(exp(aimm*(tmelt-t)) - 1) / lamc**3
  mnuccc = dum * nnuccc * pi/6*rhow * (pgam+4)(pgam+5)(pgam+6) / lamc**3
  ```
- 🎲 `dum` is the subgrid variance factor from `var_coef(relvar, 2, ...)` when `microp_uniform` is false, otherwise 1.

## 🔗🧊 Call site and limiter (`micro_pumas_v1.F90` ~line 2310)

- 📞 Inside `if (do_cldice)`: calls `immersion_freezing(...)`, then caps `nnuccc*lcldm` at the deposition ice-nuclei number `nnuccd` to prevent runaway droplet freezing (mass scaled by the same ratio), then calls `contact_freezing(...)`.

## 🧪🔀 Classical nucleation alternative (`use_hetfrz_classnuc = .true.`, ~line 2345)

- 🚫 PUMAS does not compute immersion freezing itself. It takes the input `frzimm` (1/cm3, from CAM's `hetfrz_classnuc`) and converts it:
  ```
  nnuccc = frzimm * 1.0e6 / rho
  ```
- ⚖️ Mass frozen is the mean droplet mass with limiters (number >= 1/cm3, minimum mass).

## 📊🗂️ Related

- 📈 Diagnostic `nnuccctot` in `micro_pumas_diags.F90` (~line 69).
- 🧾 CCPP metadata for the tendency in `micro_pumas_ccpp.F90` (~line 346) and `micro_pumas_ccpp.meta`.
- 🕰️ Legacy MG1 version in `micro_mg1_0.F90` (~line 1905).

## 📂 Files

/cluster/work/climate/dvillanueva/work/CESM/PUMAS/micro_pumas_utils.F90
/cluster/work/climate/dvillanueva/work/CESM/PUMAS/micro_pumas_v1.F90
