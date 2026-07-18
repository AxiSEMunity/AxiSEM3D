#!/usr/bin/env python3
"""Validate that a generated Exodus mesh has strictly positive fluid density."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def decode_elem_var_names(mesh: h5py.File) -> list[str]:
    raw = mesh["name_elem_var"][:]
    return [
        b"".join(raw[i]).decode("ascii", errors="ignore").strip().strip("\x00").strip()
        for i in range(raw.shape[0])
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", help="Path to the Exodus .e mesh file")
    args = parser.parse_args()

    path = Path(args.mesh).resolve()
    with h5py.File(path, "r") as mesh:
        names = decode_elem_var_names(mesh)
        idx_fluid = names.index("fluid") + 1
        fluid = mesh[f"vals_elem_var{idx_fluid}eb1"][0] > 0.5

        bad_counts = []
        for inode in range(4):
            idx_rho = names.index(f"RHO_{inode}") + 1
            rho = mesh[f"vals_elem_var{idx_rho}eb1"][0]
            bad = np.count_nonzero(fluid & (rho <= 0.0))
            bad_counts.append(bad)
            print(
                f"RHO_{inode}: min={rho[fluid].min():.12g}, "
                f"max={rho[fluid].max():.12g}, nonpositive_fluid_values={bad}"
            )

    total_bad = sum(bad_counts)
    if total_bad:
        raise SystemExit(
            f"Mesh validation failed: found {total_bad} nonpositive fluid density "
            f"corner values in {path}"
        )
    print(f"Mesh validation passed: all fluid RHO_* values are strictly positive in {path}")


if __name__ == "__main__":
    main()
