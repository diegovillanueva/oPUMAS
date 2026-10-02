#!/usr/bin/env python3
"""⏱️📈 Time-pressure plots of a column_loop.jl run: temperature change, cloud liquid,
cloud ice, rain and snow, plus surface precipitation, one row per column.

Usage: plot_column_loop.py column_loop.nc [--columns 3 17 ...] [--count 3]
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

PANELS = (
    ("t", "ΔT (K)", "RdBu_r"),
    ("qcn", "cloud liquid (g/kg)", "Blues"),
    ("qin", "cloud ice (g/kg)", "Purples"),
    ("qrn", "rain (g/kg)", "Greens"),
    ("qsn", "snow (g/kg)", "Oranges"),
)


def wettest_columns(ds: xr.Dataset, count: int) -> list[int]:
    """🌧️🔢 Columns with the largest mean surface precipitation."""
    mean = ds["prect"].mean("step").values
    return [int(c) for c in np.argsort(mean)[::-1][:count]]


def plot(ds: xr.Dataset, columns: list[int], figure: Path) -> None:
    """🖼️📈 One row of time-pressure panels per column."""
    hours = ds["step"].size * float(ds.attrs["deltat"]) / 3600.0
    time = np.linspace(0.0, hours, ds["step"].size)
    fig, axes = plt.subplots(
        len(columns),
        len(PANELS) + 1,
        figsize=(3.2 * (len(PANELS) + 1), 2.8 * len(columns)),
        squeeze=False,
        constrained_layout=True,
    )
    for row, col in zip(axes, columns):
        pressure = ds["p"].isel(col=col).values / 100.0
        for ax, (name, label, cmap) in zip(row, PANELS):
            field = ds[name].isel(col=col).transpose("lev", "step").values
            if name == "t":
                field = field - field[:, :1]
                limit = max(np.abs(field).max(), 1e-3)
                mesh = ax.pcolormesh(time, pressure, field, cmap=cmap, vmin=-limit, vmax=limit)
            else:
                mesh = ax.pcolormesh(time, pressure, field * 1e3, cmap=cmap)
            fig.colorbar(mesh, ax=ax, label=label)
            ax.set_ylim(pressure.max(), max(pressure.min(), 100.0))
            ax.set_xlabel("hours")
        row[0].set_ylabel(f"column {col}\npressure (hPa)")
        ax = row[-1]
        ax.plot(time, ds["prect"].isel(col=col).values * 3.6e6, label="total")
        ax.plot(time, ds["preci"].isel(col=col).values * 3.6e6, label="ice")
        ax.set_xlabel("hours")
        ax.set_ylabel("surface precipitation (mm/h)")
        ax.legend()
    fig.savefig(figure, dpi=120)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("loop", type=Path)
    parser.add_argument("--columns", type=int, nargs="*")
    parser.add_argument("--count", type=int, default=3)
    args = parser.parse_args()
    ds = xr.open_dataset(args.loop)
    columns = args.columns or wettest_columns(ds, args.count)
    figure = args.loop.with_suffix(".png").resolve()
    plot(ds, columns, figure)
    hours = ds["step"].size * float(ds.attrs["deltat"]) / 3600.0
    notes = figure.with_suffix(".md")
    notes.write_text(
        f"🖼️⏱️ Offline single-column run of the oPUMAS Julia twin: {ds['step'].size} steps "
        f"of {float(ds.attrs['deltat']):.0f} s ({hours:.0f} h), columns {columns} of "
        f"`{ds.attrs['inputs']}`, settings `{ds.attrs['settings']}`.\n"
        "🌍🔒 Host forcing (cloud fractions, activated nuclei, pressure) is held fixed; "
        "the state advances by forward Euler with the microphysics tendencies only, so "
        "this is a microphysics box model, not a forecast.\n"
        "📈 Panels: temperature change since the first step, cloud liquid, cloud ice, "
        "rain and snow mass (g/kg), surface precipitation (mm/h).\n"
    )
    print(figure)
    print(notes)


if __name__ == "__main__":
    main()
