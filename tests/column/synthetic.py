#!/usr/bin/env python3
"""🧪🌦️ Hand-made columns for processes that real CAM columns rarely trigger.

🧱 Pressure, dust and subgrid fields come from the first dumped column; temperature,
   humidity and hydrometeors are set by each recipe below.
Usage: synthetic.py <selected inputs.nc> <synthetic_inputs.nc>
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

T0 = 273.15
HYDRO = ["qcn", "qin", "ncn", "nin", "qrn", "qsn", "nrn", "nsn", "qgr", "ngr"]


def qsat(t: np.ndarray, p: np.ndarray, ice: bool = False) -> np.ndarray:
    """💧📏 Bolton-style saturation mixing ratio; only used to set relative humidity."""
    tc = t - T0
    es = 611.2 * (np.exp(22.46 * tc / (t - 0.53)) if ice else np.exp(17.67 * tc / (t - 29.65)))
    es = np.minimum(es, 0.5 * p)
    return 0.622 * es / (p - 0.378 * es)


def base(tmpl: xr.Dataset, t_sfc: float, lapse: float = 6.5e-3) -> dict[str, np.ndarray]:
    """🌡️🧊 A cloud-free column with a standard lapse rate and 50 % relative humidity."""
    c = {v: tmpl[v].values.copy() for v in tmpl.data_vars}
    p = c["p"]
    z = 7400.0 * np.log(p.max() / p)
    c["t"] = np.maximum(t_sfc - lapse * z, 200.0)
    c["q"] = 0.5 * qsat(c["t"], p)
    for v in HYDRO + ["cldn", "liqcldf", "icecldf", "tnd_qsnow", "tnd_nsnow",
                      "frzimm", "frzcnt", "frzdep", "npccn", "naai"]:
        c[v] = np.zeros_like(p)
    c["re_ice"] = np.full_like(p, 25e-6)
    return c


def layer(c: dict, lo: float, hi: float) -> np.ndarray:
    """📐🎚️ Levels with pressure between lo and hi (Pa)."""
    return (c["p"] >= lo) & (c["p"] <= hi)


def cloud(c: dict, k: np.ndarray, qc: float = 0.0, qi: float = 0.0, frac: float = 0.8) -> None:
    """☁️💧 Saturated cloud layer with liquid and/or ice and matching numbers."""
    c["cldn"][k] = frac
    c["liqcldf"][k] = frac if qc > 0 else 0.0
    c["icecldf"][k] = frac if qi > 0 else 0.0
    c["q"][k] = qsat(c["t"][k], c["p"][k], ice=c["t"][k].mean() < T0 and qc == 0)
    c["qcn"][k], c["ncn"][k] = qc, 100e6 * (qc > 0)
    c["qin"][k], c["nin"][k] = qi, 1e5 * (qi > 0)
    c["npccn"][k] = 1e5 * (qc > 0)
    c["naai"][k] = 2e5 * (qi > 0 or c["t"][k].mean() < T0)


def recipes(tmpl: xr.Dataset) -> dict[str, dict]:
    out = {}
    c = base(tmpl, 288.0)
    out["clear_dry"] = c

    c = base(tmpl, 275.0)
    k = layer(c, 60000, 80000)
    cloud(c, k, qc=3e-4, qi=0.0)
    c["qrn"][k | layer(c, 80000, 2e5)], c["nrn"][k | layer(c, 80000, 2e5)] = 2e-4, 5e3
    c["t"][k] -= 8.0  # 🧊💧 supercooled liquid and rain near -20 C: freezing, Bergeron, riming
    c["frzimm"][k], c["frzcnt"][k], c["frzdep"][k] = 1e-3, 1e-4, 1e-4
    out["supercooled_rain"] = c

    c = base(tmpl, 280.0)
    k = layer(c, 20000, 35000)
    cloud(c, k, qc=2e-4)  # ❄️🔥 liquid below -40 C: homogeneous freezing
    out["homogeneous_freeze"] = c

    c = base(tmpl, 285.0)
    k = layer(c, 40000, 75000)
    cloud(c, k, qc=1e-4, qi=1e-4)
    falling = layer(c, 40000, 2e5)
    c["qsn"][falling], c["nsn"][falling] = 5e-4, 1e4
    c["qgr"][falling], c["ngr"][falling] = 5e-4, 1e3
    c["qrn"][falling], c["nrn"][falling] = 1e-4, 1e3
    out["melting_snow_graupel"] = c

    c = base(tmpl, 300.0)
    c["q"] *= 0.4
    falling = layer(c, 50000, 2e5)
    c["qrn"][falling], c["nrn"][falling] = 3e-4, 1e4
    c["qsn"][layer(c, 50000, 60000)], c["nsn"][layer(c, 50000, 60000)] = 2e-4, 1e4
    out["dry_evaporation"] = c

    c = base(tmpl, 270.0)
    k = layer(c, 50000, 85000)
    cloud(c, k, qc=8e-4, qi=2e-4, frac=1.0)
    c["qsn"][k], c["nsn"][k] = 1e-3, 2e4
    c["tnd_qsnow"][k], c["tnd_nsnow"][k] = 1e-7, 1e-2
    out["riming_graupel"] = c

    c = base(tmpl, 298.0)
    k = layer(c, 85000, 95000)
    cloud(c, k, qc=1.5e-3, frac=1.0)
    c["ncn"][k] = 20e6  # 🌧️💧 few large droplets: strong autoconversion
    out["heavy_drizzle"] = c

    c = base(tmpl, 260.0)
    k = layer(c, 25000, 40000)
    c["q"][k] = 1.3 * qsat(c["t"][k], c["p"][k], ice=True)
    c["cldn"][k], c["icecldf"][k] = 0.3, 0.3
    c["naai"][k] = 5e5  # 🧊🌫️ ice-supersaturated, cloud-free: nucleation and deposition
    out["ice_supersat"] = c

    c = base(tmpl, 290.0)
    k = layer(c, 80000, 95000)
    cloud(c, k, qi=1e-4, frac=0.5)  # 🔥🧊 cloud ice below the freezing level: instant melting
    c["q"][k] = qsat(c["t"][k], c["p"][k])
    out["warm_ice_melt"] = c
    return out


def main() -> None:
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    ds = xr.open_dataset(src)
    # 🔄📐 Files are written by Fortran, so xarray sees (lev, col); work column-first here.
    order = {v: ds[v].dims for v in ds.data_vars}
    tmpl = ds.isel(col=[0]).transpose("col", ...)
    made = recipes(tmpl.isel(col=0))
    names = list(made)
    data = {}
    for v in tmpl.data_vars:
        rows = [np.asarray(made[n].get(v, tmpl[v].values[0])) for n in names]
        data[v] = (tmpl[v].dims, np.stack(rows))
    for v in ("source_nstep", "source_lchnk", "source_icol"):
        if v in data:
            data[v] = ("col", np.full(len(names), -1, dtype="i4"))
    data["synthetic_name"] = ("col", np.array(names, dtype=object))
    out = xr.Dataset(data)
    out = out.assign({v: out[v].transpose(*order[v]) for v in order})
    out.attrs = ds.attrs
    out.to_netcdf(dst)
    print(dst.resolve())


if __name__ == "__main__":
    main()
