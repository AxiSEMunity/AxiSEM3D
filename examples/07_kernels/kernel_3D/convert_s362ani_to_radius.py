"""
convert_s362ani_to_radius.py
============================
Convert S362ANI_percent.nc from depth-based to radius-based coordinates,
as required by AxiSEM3D's StructuredGridV3D class (vertical axis must be
monotonically *increasing* radius).  Also converts perturbation values
from percent to fractional form (divide by 100), because AxiSEM3D's
REF1D mode applies: v_3d = v_1d * (1 + factor * data).

Source: S362ANI_percent.nc  (depth axis in km, shallow → deep, percent)
Output: input_forward/S362ANI_radius.nc  (radius axis in m, deep → shallow, fractional)

Run::

    python convert_s362ani_to_radius.py
"""

import os
import numpy as np
import netCDF4 as nc

# ── paths ────────────────────────────────────────────────────────────────────
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))

# Resolve source relative to the repository root
_REPO_ROOT   = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", ".."))
SOURCE_FILE  = os.path.join(
    _REPO_ROOT, "AxiSEM3D_Kernels",
    "examples", "data", "3D_KERNEL_EXAMPLE_30s", "input", "S362ANI_percent.nc",
)
OUTPUT_FILE  = os.path.join(_SCRIPT_DIR, "input_forward", "S362ANI_radius.nc")

# Earth radius used for conversion (m)
EARTH_RADIUS_M = 6_371_000.0


def convert_s362ani(source_file=SOURCE_FILE, output_file=OUTPUT_FILE):
    """
    Read S362ANI_percent.nc, convert depth→radius, flip axis to be
    monotonically increasing, and write S362ANI_radius.nc.

    Parameters
    ----------
    source_file : str
        Path to the original depth-based NetCDF file.
    output_file : str
        Destination path for radius-based output.

    Returns
    -------
    str
        Path to the written output file.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with nc.Dataset(source_file, "r") as src:
        depth_km  = src["depth"][:]          # shape (25,)  km, ascending
        latitude  = src["latitude"][:]       # shape (91,)  degrees
        longitude = src["longitude"][:]      # shape (181,) degrees
        dvs       = src["dvs"][:]            # shape (25, 91, 181)
        dvsv      = src["dvsv"][:]
        dvsh      = src["dvsh"][:]

    # depth (km, ascending) → radius (m): radius = R_E - depth*1000
    # depth ascending → radius descending → flip to make radius ascending
    radius_m = EARTH_RADIUS_M - depth_km * 1000.0   # still in original order
    radius_m = radius_m[::-1].copy()                # flip: now ascending

    # Flip data variables along axis 0 to match
    dvs  = dvs [::-1].copy()
    dvsv = dvsv[::-1].copy()
    dvsh = dvsh[::-1].copy()

    # Convert percent → fractional (AxiSEM3D REF1D: v = v_1d * (1 + factor * data))
    dvs  /= 100.0
    dvsv /= 100.0
    dvsh /= 100.0

    with nc.Dataset(output_file, "w", format="NETCDF4") as ds:
        ds.description = (
            "S362ANI radial anisotropy model, converted from depth (km) to "
            "radius (m) for use with AxiSEM3D StructuredGridV3D. "
            "Perturbations are fractional relative to PREM (divided by 100)."
        )
        ds.source = f"Converted from {os.path.basename(source_file)}"

        # Dimensions
        ds.createDimension("radius",    len(radius_m))
        ds.createDimension("latitude",  len(latitude))
        ds.createDimension("longitude", len(longitude))

        # Coordinate variables
        r_var = ds.createVariable("radius",    "f8", ("radius",))
        r_var.units     = "m"
        r_var.long_name = "Radius from Earth centre"
        r_var[:]        = radius_m

        lat_var = ds.createVariable("latitude",  "f8", ("latitude",))
        lat_var.units = "degrees_north"
        lat_var[:]    = latitude

        lon_var = ds.createVariable("longitude", "f8", ("longitude",))
        lon_var.units = "degrees_east"
        lon_var[:]    = longitude

        # Data variables  (radius, latitude, longitude)
        for name, data in [("dvs", dvs), ("dvsv", dvsv), ("dvsh", dvsh)]:
            v = ds.createVariable(name, "f8", ("radius", "latitude", "longitude"))
            v.units     = "fractional"
            v.long_name = f"Relative S-velocity perturbation ({name}), fractional"
            v[:]        = data

    print(f"Written: {output_file}")
    print(f"  radius : {radius_m[0]:.0f} m .. {radius_m[-1]:.0f} m  "
          f"({len(radius_m)} pts) — monotonically "
          f"{'increasing' if float(radius_m[-1]) > float(radius_m[0]) else 'DECREASING'}")
    print(f"  lat    : {float(latitude[0])}° .. {float(latitude[-1])}°")
    print(f"  lon    : {float(longitude[0])}° .. {float(longitude[-1])}°")
    return output_file


if __name__ == "__main__":
    convert_s362ani()
