#!/usr/bin/env python3
"""
Example 06: Finite Source Local (San Francisco Bay Area) — publication figures.
  Fig 1: Velocity wavefield snapshots (4-panel)
  Fig 2: Peak ground velocity map
  Fig 3: Velocity seismograms at high-PGV stations
"""
import os
import glob
import netCDF4
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from pub_style import apply_style, save_fig, COLORS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature

apply_style()
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['pdf.use14corefonts'] = True    # reference core Helvetica (Illustrator-editable)
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['axes.unicode_minus'] = False

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
NC_DIR = os.path.join(DATA_DIR, 'stations')

if not os.path.isdir(NC_DIR):
    print('Ex06: Data not yet available.')
    import sys; sys.exit(0)


def read_stations():
    """Combine the per-rank NETCDF station files into one set of arrays.
    dask is unavailable, so open each rank file with netCDF4 and concat the
    station dimension manually (no xr.open_mfdataset)."""
    fns = sorted(glob.glob(os.path.join(NC_DIR, 'axisem3d_synthetics.nc.*')))
    names, disp_list, time = [], [], None
    for f in fns:
        ds = netCDF4.Dataset(f)
        ls = ds.variables['list_station'][:]
        for srow in ls:
            s = b''.join(srow).decode('utf-8', 'ignore').strip().strip('\x00')
            names.append(s)
        if time is None:
            time = np.array(ds.variables['data_time'][:])
        disp_list.append(np.array(ds.variables['data_wave'][:]))
        ds.close()
    return np.array(names), time, np.concatenate(disp_list, axis=0)


station_names_decoded, t_coords, all_displacement = read_stations()

stations = pd.read_csv(
    os.path.join(DATA_DIR, 'STATIONS_OUTPUT.txt'),
    sep=r'\s+', header=2,
    names=['name', 'network', 'x', 'y', 'useless', 'depth'])

stations['name_network'] = [
    f'{net}.{nam:04d}' for nam, net in zip(stations['name'], stations['network'])]

permutation_vector = np.array(
    [(nn == station_names_decoded).argmax() for nn in stations['name_network']])
stations['permutation'] = permutation_vector

lat_pole, lon_pole = 37.7, -122.1
x_utm_origin, y_utm_origin = ccrs.UTM(zone=11).transform_point(
    lon_pole, lat_pole, src_crs=ccrs.Geodetic())
stations['x_utm'] = stations['x'] + x_utm_origin
stations['y_utm'] = stations['y'] + y_utm_origin

x_coords = stations['x_utm'].values
y_coords = stations['y_utm'].values
displacement_data = all_displacement[stations['permutation'], :, :]

dt = t_coords[1] - t_coords[0]
velocity_data = np.gradient(displacement_data, dt, axis=2)

velocity_norm = np.sqrt(np.sum(velocity_data**2, axis=1))
pgv_final = np.max(velocity_norm, axis=1)

# regular grid for continuous-field interpolation, spanning the station extent;
# griddata(method='linear') returns NaN outside the stations' convex hull, so
# the field is only ever drawn where the simulation actually sampled it.
GRID_N = 300
gx = np.linspace(x_coords.min(), x_coords.max(), GRID_N)
gy = np.linspace(y_coords.min(), y_coords.max(), GRID_N)
grid_x, grid_y = np.meshgrid(gx, gy)


def interp_field(values):
    return griddata((x_coords, y_coords), values, (grid_x, grid_y), method='linear')

# ==== Figure 1: Wavefield snapshots ====
print('Ex06 Fig 1: Velocity wavefield snapshots')
t_frac = 0.6
tstep = min(int(t_frac * len(t_coords)), len(t_coords) - 1)

proj = ccrs.UTM(zone=11)
vel_max = float(np.percentile(np.abs(velocity_data), 99.5))
# PGV (a running/final MAX over time) is systematically larger than any single
# instantaneous-velocity snapshot, so sharing vel_max saturates most of the
# PGV panel; give PGV its own scale from its own distribution.
pgv_max = float(np.percentile(pgv_final, 99.0))

fig, axes = plt.subplots(2, 2, figsize=(8, 8.8), dpi=150,
                          subplot_kw={'projection': proj})
titles = ['Vertical (Z)', 'East (E)', 'North (N)', 'PGV to time']
channels = [2, 0, 1, None]
cmaps = ['RdBu_r', 'RdBu_r', 'RdBu_r', 'YlOrRd']
pgv_running = np.max(velocity_norm[:, :tstep + 1], axis=1)

for idx, ax in enumerate(axes.flat):
    ax.set_extent((-122.7, -121.5, 37.2, 38.2), crs=ccrs.Geodetic())
    try:
        ax.coastlines(resolution='10m', linewidth=0.9, color='0.25')
    except Exception:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.9, edgecolor='0.25')
    ax.set_title(titles[idx], fontsize=10)

    if idx < 3:
        field = interp_field(velocity_data[:, channels[idx], tstep])
        sc = ax.pcolormesh(grid_x, grid_y, field, vmin=-vel_max, vmax=vel_max,
                           cmap=cmaps[idx], shading='auto', rasterized=True,
                           transform=proj)
    else:
        field = interp_field(pgv_running)
        sc = ax.pcolormesh(grid_x, grid_y, field, vmin=0, vmax=pgv_max,
                           cmap=cmaps[idx], shading='auto', rasterized=True,
                           transform=proj)

cbar_vel = fig.colorbar(
    plt.cm.ScalarMappable(norm=plt.Normalize(-vel_max, vel_max), cmap='RdBu_r'),
    ax=axes[0, :].tolist(), shrink=0.6, pad=0.02, orientation='horizontal',
    location='top')
cbar_vel.set_label('Velocity (m/s)', fontsize=9)

cbar_pgv = fig.colorbar(
    plt.cm.ScalarMappable(norm=plt.Normalize(0, pgv_max), cmap='YlOrRd'),
    ax=axes[1, 1], shrink=0.7, pad=0.02, orientation='horizontal',
    location='bottom')
cbar_pgv.set_label('PGV (m/s)', fontsize=9)

fig.suptitle(f'SFBA Finite Fault — t = {t_coords[tstep]:.1f} s', fontsize=12)
save_fig(fig, os.path.join(SCRIPT_DIR, 'ex06_wavefield_snapshots'))

# ==== Figure 2: PGV map ====
print('Ex06 Fig 2: Peak ground velocity map')
fig, ax = plt.subplots(figsize=(6, 5), dpi=150,
                        subplot_kw={'projection': proj})
ax.set_extent((-122.7, -121.5, 37.2, 38.2), crs=ccrs.Geodetic())
try:
    ax.coastlines(resolution='10m', linewidth=1.0, color='0.2')
except Exception:
    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='0.2')

pgv_field = interp_field(pgv_final)
sc = ax.pcolormesh(grid_x, grid_y, pgv_field, vmin=0, vmax=pgv_max,
                   cmap='YlOrRd', shading='auto', rasterized=True,
                   transform=proj)
cbar = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.03)
cbar.set_label('Peak Ground Velocity (m/s)', fontsize=10)
ax.set_title('SFBA Finite Fault — PGV', fontsize=11)
save_fig(fig, os.path.join(SCRIPT_DIR, 'ex06_pgv_map'))

# ==== Figure 3: Velocity seismograms ====
print('Ex06 Fig 3: Velocity seismograms')
n_sel = 5
np.random.seed(42)
high_pgv_idx = np.argsort(pgv_final)[-50:]
sel_indices = np.sort(np.random.choice(high_pgv_idx, n_sel, replace=False))

fig, axes = plt.subplots(n_sel, 1, figsize=(8, 2 * n_sel), sharex=True, sharey=True)
ch_labels = ['E', 'N', 'Z']
ch_colors = [COLORS['primary'], COLORS['accent'], COLORS['secondary']]

for i, si in enumerate(sel_indices):
    ax = axes[i]
    for ch in range(3):
        ax.plot(t_coords, velocity_data[si, ch, :],
                color=ch_colors[ch], lw=0.7, label=ch_labels[ch], alpha=0.85)
    ax.set_ylabel('Vel. (m/s)', fontsize=8)
    sx, sy = stations.iloc[si]['x'] / 1e3, stations.iloc[si]['y'] / 1e3
    ax.text(0.02, 0.92, f'({sx:.0f}, {sy:.0f}) km',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round,pad=0.2', fc='wheat', alpha=0.7))
    if i == 0:
        ax.legend(loc='upper right', fontsize=7, ncol=3)

axes[-1].set_xlabel('Time (s)')
axes[0].set_title('SFBA Finite Fault — Velocity at high-PGV stations', fontsize=11)
plt.tight_layout()
save_fig(fig, os.path.join(SCRIPT_DIR, 'ex06_velocity_seismograms'))

print('Ex06 figures complete.')
