#!/usr/bin/env python3
"""🔎🧊 Pick about 200 columns from a CAM column dump so every process is exercised.

📥 Reads inputs.n*.c*.nc / outputs.n*.c*.nc written by the dump case (tests/column/dump).
📤 Writes inputs.nc, cam_outputs.nc and coverage.csv into the output directory.
"""

import argparse
from pathlib import Path

import numpy as np
import xarray as xr

QSMALL = 1e-18
MINCLD = 1e-4


def load(dump: Path, kind: str) -> xr.Dataset:
    """🧷📚 Concatenate all dump files of one kind along col, tagging each column."""
    parts = []
    for f in sorted(dump.glob(f"{kind}.n*.c*.nc")):
        ds = xr.open_dataset(f).load()
        n = ds.sizes["col"]
        ds["source_nstep"] = ("col", np.full(n, int(f.name.split(".")[1][1:]), dtype="i4"))
        ds["source_lchnk"] = ("col", np.full(n, int(f.name.split(".")[2][1:]), dtype="i4"))
        ds["source_icol"] = ("col", np.arange(1, n + 1, dtype="i4"))
        parts.append(ds.drop_attrs() if hasattr(ds, "drop_attrs") else ds)
    if not parts:
        raise SystemExit(f"no {kind} files in {dump}")
    out = xr.concat(parts, dim="col", data_vars="all", coords="minimal")
    out.attrs = {"deltat": parts[0].attrs.get("deltat", np.nan),
                 "nlev": parts[0].attrs.get("nlev", parts[0].sizes["lev"])}
    first = xr.open_dataset(sorted(dump.glob(f"{kind}.n*.c*.nc"))[0])
    out.attrs.update({k: v for k, v in first.attrs.items() if k not in ("nstep", "lchnk")})
    return out


def regimes(inp: xr.Dataset) -> dict[str, np.ndarray]:
    """🌦️🏷️ Boolean masks over columns for the cloud regimes named in the plan."""
    t, cld = inp["t"], inp["cldn"]
    qc, qi = inp["qcn"], inp["qin"]
    precip = inp["qrn"] + inp["qsn"] + inp["qgr"]
    liq = qc > 1e-6
    ice = qi > 1e-7
    return {
        "warm_cloud": ((liq & (t > 273.15)).any("lev") & ~(ice.any("lev"))).values,
        "mixed_phase": (liq & ice & (t < 273.15)).any("lev").values,
        "cirrus": ((ice & (t < 235.0)).any("lev") & ~(liq.any("lev"))).values,
        "deep_outflow": ((ice & (inp["p"] < 30000)).any("lev") & (precip > 1e-4).any("lev")).values,
        "polar": (np.abs(inp["lat"]) > 60).values,
        "clear_sky": (cld.max("lev") < MINCLD).values,
        "precip_no_cloud": ((precip > QSMALL) & (cld < MINCLD)).any("lev").values,
    }


def active(out: xr.Dataset) -> dict[str, np.ndarray]:
    """⚙️🔢 Per process rate: which columns have a non-zero value anywhere."""
    rates = [v for v in out.data_vars if v.startswith("proc_rates_") and "lev" in out[v].dims]
    return {v: (out[v] != 0).any("lev").values for v in rates}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dump", type=Path)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--per-process", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20261002)
    a = ap.parse_args()

    inp, out = load(a.dump, "inputs"), load(a.dump, "outputs")
    ncol = inp.sizes["col"]
    act, reg = active(out), regimes(inp)
    rng = np.random.default_rng(a.seed)
    chosen: list[int] = []

    def take(mask: np.ndarray, k: int) -> None:
        have = int(mask[chosen].sum()) if chosen else 0
        pool = [i for i in rng.permutation(np.flatnonzero(mask)) if i not in chosen]
        chosen.extend(pool[: max(k - have, 0)])

    # 🧮🎯 Rare processes first, so common ones are covered by the same columns.
    for name in sorted(act, key=lambda v: act[v].sum()):
        take(act[name], a.per_process)
    for mask in reg.values():
        take(mask, 15)
    take(np.ones(ncol, bool), a.n)
    idx = np.sort(np.array(chosen[: max(a.n, len(chosen))]))

    a.outdir.mkdir(parents=True, exist_ok=True)
    inp.isel(col=idx).to_netcdf(a.outdir / "inputs.nc")
    out.isel(col=idx).to_netcdf(a.outdir / "cam_outputs.nc")
    with open(a.outdir / "coverage.csv", "w") as f:
        f.write("name,selected_nonzero,dump_nonzero,dump_columns\n")
        for name, m in sorted(reg.items()) + sorted(act.items()):
            f.write(f"{name.removeprefix('proc_rates_')},{int(m[idx].sum())},{int(m.sum())},{ncol}\n")
    print(f"selected {idx.size} of {ncol} columns")
    for p in ("inputs.nc", "cam_outputs.nc", "coverage.csv"):
        print((a.outdir / p).resolve())


if __name__ == "__main__":
    main()
