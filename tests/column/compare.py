#!/usr/bin/env python3
"""⚖️🔢 Compare two column output files variable by variable.

🎯 Exit 0 when every shared variable agrees within --rtol (relative) or --atol (absolute),
   1 otherwise. --rtol 0 --atol 0 asks for bit-for-bit agreement.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import xarray as xr


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("new", type=Path)
    ap.add_argument("ref", type=Path)
    ap.add_argument("--rtol", type=float, default=0.0)
    ap.add_argument("--atol", type=float, default=0.0)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    new, ref = xr.open_dataset(a.new), xr.open_dataset(a.ref)
    bad, worst = [], (0.0, "", -1)
    shared = [v for v in ref.data_vars if v in new.data_vars and not v.startswith("source_")]
    for v in shared:
        x, y = new[v].values.astype("f8"), ref[v].values.astype("f8")
        both_nan = np.isnan(x) & np.isnan(y)
        diff = np.where(both_nan, 0.0, np.abs(x - y))
        rel = np.where(y != 0, diff / np.abs(np.where(y != 0, y, 1)), np.where(diff > 0, np.inf, 0))
        ok = (diff <= a.atol) | (rel <= a.rtol)
        ok &= ~(np.isnan(x) ^ np.isnan(y))
        r = float(np.max(np.where(diff > a.atol, rel, 0))) if rel.size else 0.0
        if r > worst[0]:
            col = int(np.unravel_index(np.argmax(np.where(diff > a.atol, rel, 0)), rel.shape)[0])
            worst = (r, v, col)
        if not ok.all():
            bad.append((v, int((~ok).sum()), r))
    missing = [v for v in ref.data_vars if v not in new.data_vars and not v.startswith("source_")]
    if not a.quiet:
        for v, n, r in bad:
            print(f"FAIL {v}: {n} values off, max rel {r:.3e}")
        for v in missing:
            print(f"MISSING {v}")
    print(f"{'PASS' if not bad and not missing else 'FAIL'} {len(shared)} variables; "
          f"worst rel {worst[0]:.3e} in {worst[1] or '-'} col {worst[2]}")
    sys.exit(0 if not bad and not missing else 1)


if __name__ == "__main__":
    main()
