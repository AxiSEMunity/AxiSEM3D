"""
create_moho_topography.py
=========================
Generate a Gaussian Moho undulation NetCDF file compatible with
AxiSEM3D's StructuredGridG3D model class.

Grid:   1° × 1° global (lat -90..90, lon -180..180)
Model:  Gaussian bump centred at (lat=0°, lon=20°),
        amplitude = 5000 m, angular sigma ≈ 30°.

Run::

    python create_moho_topography.py

Output: input_forward/moho_topography.nc
"""

import os
import numpy as np
import netCDF4 as nc

# ── configuration ─────────────────────────────────────────────────────────────
OUTPUT_PATH = os.path.join(os.path.dirname(__file__),
                           "input_forward", "moho_topography.nc")

# Gaussian parameters
AMPLITUDE_M  = 5000.0          # metres
CENTER_LAT   = 0.0             # degrees
CENTER_LON   = 20.0            # degrees
# σ such that the Gaussian drops to 1/10 of its peak at ~30° angular distance
# exp(-0.5 * (30/σ)²) = 0.1  →  σ = 30 / sqrt(2 * ln(10)) ≈ 13.97°
SIGMA_DEG    = 30.0 / np.sqrt(2.0 * np.log(10.0))

# NetCDF variable names expected by StructuredGridG3D
FACTOR       = 1.0             # multiplied by undulation values


def _angular_distance_deg(lat1, lon1, lat2, lon2):
    """Great-circle angular distance (degrees) using the haversine formula."""
    lat1_r, lat2_r = np.radians(lat1), np.radians(lat2)
    dlon_r = np.radians(lon2 - lon1)
    a = (np.sin((lat2_r - lat1_r) / 2.0) ** 2
         + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon_r / 2.0) ** 2)
    return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))


def create_moho_topography(output_path=OUTPUT_PATH,
                           amplitude_m=AMPLITUDE_M,
                           center_lat=CENTER_LAT,
                           center_lon=CENTER_LON,
                           sigma_deg=SIGMA_DEG,
                           factor=FACTOR):
    """
    Generate and write the Moho undulation NetCDF file.

    Parameters
    ----------
    output_path : str
        Destination file path.
    amplitude_m : float
        Peak undulation amplitude in metres.
    center_lat : float
        Latitude of Gaussian centre (degrees).
    center_lon : float
        Longitude of Gaussian centre (degrees).
    sigma_deg : float
        Angular standard deviation of the Gaussian (degrees).
    factor : float
        Multiplicative factor applied to undulation values (YAML ``factor``).
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # 1° × 1° global grid
    latitudes  = np.arange(-90.0, 91.0, 1.0)   # 181 points
    longitudes = np.arange(-180.0, 181.0, 1.0)  # 361 points

    lon2d, lat2d = np.meshgrid(longitudes, latitudes)   # shape (181, 361)

    # Angular distance from Gaussian centre (degrees)
    dist_deg = _angular_distance_deg(lat2d, lon2d, center_lat, center_lon)

    # Gaussian undulation field (metres)  — multiplied by factor
    undulation = factor * amplitude_m * np.exp(
        -0.5 * (dist_deg / sigma_deg) ** 2
    )

    # Write NetCDF
    with nc.Dataset(output_path, "w", format="NETCDF4") as ds:
        ds.description = (
            "Moho topographic undulation for AxiSEM3D StructuredGridG3D. "
            f"Gaussian bump: amplitude={amplitude_m} m, "
            f"centre=({center_lat}°N, {center_lon}°E), "
            f"sigma={sigma_deg:.4f}°."
        )

        # Dimensions
        ds.createDimension("latitude",  len(latitudes))
        ds.createDimension("longitude", len(longitudes))

        # Coordinate variables
        lat_var = ds.createVariable("latitude",  "f8", ("latitude",))
        lat_var.units = "degrees_north"
        lat_var[:] = latitudes

        lon_var = ds.createVariable("longitude", "f8", ("longitude",))
        lon_var.units = "degrees_east"
        lon_var[:] = longitudes

        # Undulation data variable
        und_var = ds.createVariable("undulation_MOHO", "f8",
                                    ("latitude", "longitude"))
        und_var.units = "m"
        und_var.long_name = "Moho undulation (positive = upward)"
        und_var[:] = undulation

    print(f"Written: {output_path}")
    print(f"  lat  : {latitudes[0]}° .. {latitudes[-1]}°  ({len(latitudes)} pts)")
    print(f"  lon  : {longitudes[0]}° .. {longitudes[-1]}°  ({len(longitudes)} pts)")
    print(f"  max undulation : {undulation.max():.1f} m")
    print(f"  min undulation : {undulation.min():.1f} m")
    return output_path


if __name__ == "__main__":
    create_moho_topography()
