#!/usr/bin/env python3
"""🧪♻️ Compare two CAM restart files variable by variable, bit for bit.

🧭📌 `cmp` only works within one case, because the header holds the case name.
🔎📄 Here every variable is compared as raw bytes; `hrestpath` and global attributes are ignored.
"""

import argparse
import sys
from pathlib import Path

import netCDF4
import numpy as np

IGNORED = {"hrestpath"}  # holds the path of the history restart file, which contains the case name


def differing_variables(path_a: Path, path_b: Path) -> tuple[int, list[str]]:
    """Return the number of compared variables and the names of those that differ."""
    different: list[str] = []
    with netCDF4.Dataset(path_a) as a, netCDF4.Dataset(path_b) as b:
        a.set_auto_maskandscale(False)
        b.set_auto_maskandscale(False)
        names = sorted(set(a.variables) | set(b.variables))
        for name in names:
            if name in IGNORED:
                continue
            if name not in a.variables or name not in b.variables:
                different.append(f"{name} (only in one file)")
                continue
            va, vb = np.asarray(a[name][...]), np.asarray(b[name][...])
            if va.shape != vb.shape or va.dtype != vb.dtype or va.tobytes() != vb.tobytes():
                different.append(name)
    return len(names) - len(IGNORED & set(names)), different


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file_a", type=Path)
    parser.add_argument("file_b", type=Path)
    args = parser.parse_args()
    n, different = differing_variables(args.file_a, args.file_b)
    if different:
        print(f"DIFFERENT: {len(different)} of {n} variables differ: {', '.join(different[:12])}")
        return 1
    print(f"IDENTICAL: {n} variables bit-for-bit")
    return 0


if __name__ == "__main__":
    sys.exit(main())
