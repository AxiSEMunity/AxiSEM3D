#!/usr/bin/env python3
"""Prepare a mesher-safe .bm file for atmosphere-coupled examples.

This rewrites the atmosphere branch of a SalvusMeshLite background model onto
a denser radius grid and interpolates density in log-space so the profile
stays strictly positive. The goal is to prevent nonphysical negative density
values from being introduced into the generated Exodus mesh.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np


def split_bm(path: Path):
    header_lines = []
    data_lines = []
    with path.open("r", encoding="utf-8-sig") as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                header_lines.append(line)
                continue
            token = stripped.split()[0]
            try:
                float(token)
            except ValueError:
                header_lines.append(line)
            else:
                data_lines.append(line)
    if not data_lines:
        raise RuntimeError(f"No numeric rows found in {path}")
    return header_lines, data_lines


def parse_columns(header_lines: list[str]) -> list[str]:
    for line in header_lines:
        stripped = line.strip()
        if stripped.upper().startswith("COLUMNS "):
            return stripped.split()[1:]
    raise RuntimeError("Could not find a COLUMNS header in the .bm file.")


def parse_rows(data_lines: list[str]) -> np.ndarray:
    rows = [[float(tok) for tok in line.split()] for line in data_lines]
    return np.array(rows, dtype=float)


def get_shear_column_indices(columns: list[str]) -> list[int]:
    shear_names = [name for name in ("vs", "vsv", "vsh") if name in columns]
    if not shear_names:
        raise RuntimeError("Could not identify a shear-velocity column in the .bm file.")
    return [columns.index(name) for name in shear_names]


def find_atmosphere_start(rows: np.ndarray, shear_cols: list[int]) -> int:
    zero_shear = np.all(np.isclose(rows[:, shear_cols], 0.0), axis=1)
    tail_start = None
    for idx in range(len(rows) - 1, -1, -1):
        if zero_shear[idx]:
            tail_start = idx
        else:
            break
    if tail_start is None:
        raise RuntimeError("Could not find a zero-shear atmosphere branch.")
    return tail_start


def interpolate_segment(left: np.ndarray,
                        right: np.ndarray,
                        radius_idx: int,
                        rho_idx: int,
                        max_step: float) -> list[np.ndarray]:
    r0 = left[radius_idx]
    r1 = right[radius_idx]
    if math.isclose(r0, r1):
        return [right.copy()]

    span = r1 - r0
    nseg = max(int(math.ceil(abs(span) / max_step)), 1)
    radii = np.linspace(r0, r1, nseg + 1)[1:]

    out = []
    positive_left = left[rho_idx] > 0.0
    positive_right = right[rho_idx] > 0.0
    for radius in radii:
        t = (radius - r0) / span
        row = left + (right - left) * t
        row[radius_idx] = radius
        if positive_left and positive_right:
            row[rho_idx] = math.exp(
                math.log(left[rho_idx]) * (1.0 - t) + math.log(right[rho_idx]) * t
            )
        if row[rho_idx] <= 0.0:
            raise RuntimeError(
                f"Interpolated nonpositive density at radius {radius:.6f} m. "
                "Decrease --max-atmosphere-step or inspect the source profile."
            )
        out.append(row)
    return out


def densify_atmosphere(rows: np.ndarray,
                       columns: list[str],
                       max_step: float) -> np.ndarray:
    radius_idx = columns.index("radius")
    rho_idx = columns.index("rho")
    shear_cols = get_shear_column_indices(columns)
    atm_start = find_atmosphere_start(rows, shear_cols)

    solid_rows = [row.copy() for row in rows[:atm_start]]
    atmosphere_rows = rows[atm_start:]
    if np.any(atmosphere_rows[:, rho_idx] <= 0.0):
        raise RuntimeError("The source atmosphere branch must be strictly positive in density.")

    dense_rows = [atmosphere_rows[0].copy()]
    for left, right in zip(atmosphere_rows[:-1], atmosphere_rows[1:]):
        dense_rows.extend(interpolate_segment(left, right, radius_idx, rho_idx, max_step))

    return np.vstack(solid_rows + dense_rows)


def format_row(row: np.ndarray) -> str:
    return "    " + "  ".join(f"{value:.12g}" for value in row) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Input .bm file")
    parser.add_argument("--output", required=True, help="Output .bm file")
    parser.add_argument(
        "--max-atmosphere-step",
        type=float,
        default=100.0,
        help="Maximum radial spacing in meters used when densifying the atmosphere branch",
    )
    args = parser.parse_args()

    in_path = Path(args.input).resolve()
    out_path = Path(args.output).resolve()

    header_lines, data_lines = split_bm(in_path)
    columns = parse_columns(header_lines)
    rows = parse_rows(data_lines)
    dense_rows = densify_atmosphere(rows, columns, args.max_atmosphere_step)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for line in header_lines:
            f.write(line)
        for row in dense_rows:
            f.write(format_row(row))

    rho_idx = columns.index("rho")
    print(f"Wrote {out_path}")
    print(f"Original rows: {rows.shape[0]}")
    print(f"Prepared rows: {dense_rows.shape[0]}")
    print(
        f"Atmosphere density range: [{dense_rows[:, rho_idx].min():.12g}, "
        f"{dense_rows[:, rho_idx].max():.12g}]"
    )


if __name__ == "__main__":
    main()
