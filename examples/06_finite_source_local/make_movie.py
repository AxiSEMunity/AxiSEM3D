#!/usr/bin/env python3
"""
Example 06 SFBA finite fault -- vertical-velocity wavefield animation.

Reuses the same station read + continuous-field interpolation as make_figures.py
(scipy.griddata, masked to the stations' convex hull -- only ever draws the field
where the simulation actually sampled it). Animates the full record. Writes
../../movies/ex06_sfba_velocity_Z.mp4 (ffmpeg-free, via imageio_ffmpeg).
"""
import os
import glob
import netCDF4
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from pub_style import apply_style, COLORS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature

apply_style()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')
NC_DIR = os.path.join(DATA_DIR, 'stations')

if not os.path.isdir(NC_DIR):
    print('Ex06: Data not yet available.')
    import sys; sys.exit(0)


def read_stations():
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
proj = ccrs.UTM(zone=11)
x_utm_origin, y_utm_origin = proj.transform_point(lon_pole, lat_pole, src_crs=ccrs.Geodetic())
stations['x_utm'] = stations['x'] + x_utm_origin
stations['y_utm'] = stations['y'] + y_utm_origin

x_coords = stations['x_utm'].values
y_coords = stations['y_utm'].values
displacement_data = all_displacement[stations['permutation'], :, :]
dt = t_coords[1] - t_coords[0]
velocity_z = np.gradient(displacement_data[:, 2, :], dt, axis=1)

GRID_N = 160
gx = np.linspace(x_coords.min(), x_coords.max(), GRID_N)
gy = np.linspace(y_coords.min(), y_coords.max(), GRID_N)
grid_x, grid_y = np.meshgrid(gx, gy)
vel_max = float(np.percentile(np.abs(velocity_z), 99.5))

fig = plt.figure(figsize=(6, 6.4))
ax = fig.add_subplot(projection=proj)
ax.set_extent((-122.7, -121.5, 37.2, 38.2), crs=ccrs.Geodetic())
try:
    ax.coastlines(resolution='10m', linewidth=1.0, color='0.25')
except Exception:
    ax.add_feature(cfeature.COASTLINE, linewidth=1.0, edgecolor='0.25')
ttl = ax.set_title('', fontsize=13)
state = {'mesh': None}


def update(fr):
    field = griddata((x_coords, y_coords), velocity_z[:, fr], (grid_x, grid_y), method='linear')
    if state['mesh'] is not None:
        state['mesh'].remove()
    state['mesh'] = ax.pcolormesh(grid_x, grid_y, field, vmin=-vel_max, vmax=vel_max,
                                  cmap='RdBu_r', shading='auto', rasterized=True,
                                  transform=proj, zorder=1)
    ttl.set_text(f'SFBA Finite Fault — vertical velocity   t = {t_coords[fr]:.1f} s')
    return state['mesh'], ttl


frames = range(len(t_coords))
print(f'{len(frames)} frames, t = {t_coords[0]:.1f}..{t_coords[-1]:.1f} s')

import imageio_ffmpeg
plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
ani = animation.FuncAnimation(fig, update, frames=frames, blit=False)
mdir = os.path.join(SCRIPT_DIR, '..', '..', 'movies')
os.makedirs(mdir, exist_ok=True)
out = os.path.join(mdir, 'ex06_sfba_velocity_Z.mp4')
ani.save(out, writer=animation.FFMpegWriter(fps=12, bitrate=6000), dpi=120)
print('Saved', out)
