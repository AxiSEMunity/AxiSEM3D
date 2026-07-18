#!/usr/bin/env python3
"""
Example 00: Global 1D (PREM) -- publication figures.

  Fig 1: Station map with the event (Cartopy, coastlines only)
  Fig 2: Three-component (R/T/Z) displacement seismograms at IU.ANMO
  Fig 3: USArray wavefield time-lapse grids (U3 vertical displacement
         and R3 vertical rotation)

This script reads the raw AxiSEM3D output that THIS example produces under
./output (set by INPUT below). Run the simulation first:

    mpirun -np 4 axisem3d --input input --output output

Event location and station coordinates are read from ./input (committed with
the example). Figures are written to ./figures as both PDF and PNG.
"""
import sys
import os
import numpy as np
import yaml

# ---- locate this example and the shared style helper --------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
# pub_style.py lives one level up, at the examples root.
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import (apply_style, save_fig, epicentral_distance_deg,
                       wavefield_timelapse_map, CHANNEL_COLORS, CHANNEL_LABELS,
                       COLORS)
import matplotlib.pyplot as plt

apply_style()

# ---- input / output locations -------------------------------------------- #
# INPUT: the AxiSEM3D run output this example produces. Override if you wrote
# the solver output somewhere else (e.g. INPUT = '/path/to/output').
INPUT = os.path.join(HERE, 'output')
INPUT_DIR = os.path.join(HERE, 'input')        # event + station coordinates
FIG_DIR = os.path.join(HERE, 'figures')

# ---- read event ---------------------------------------------------------- #
with open(os.path.join(INPUT_DIR, 'inparam.source.yaml')) as f:
    src = yaml.safe_load(f)
loc = src['list_of_sources'][0]['VIRGINIA_201108231751A']['location']
ev_lat, ev_lon = loc['latitude_longitude']

# ---- read GSN station info ----------------------------------------------- #
info_GSN = np.loadtxt(os.path.join(INPUT_DIR, 'GSN.txt'), dtype=str, skiprows=3)
st_lat = info_GSN[:, 2].astype(float)
st_lon = info_GSN[:, 3].astype(float)

# ==== Figure 1: Station map ============================================== #
print('Ex00 Fig 1: Station map')
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    fig = plt.figure(figsize=(7, 3.5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=-77))
    ax.set_global()
    # White background, coarse grey coastlines only (no land/ocean fill).
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4,
                   color=COLORS['gray'])
    ax.scatter(st_lon, st_lat, transform=ccrs.PlateCarree(),
               s=18, c=COLORS['primary'], marker='^', linewidths=0.3,
               edgecolors='k', zorder=10, label='GSN stations')
    ax.scatter(ev_lon, ev_lat, transform=ccrs.PlateCarree(),
               s=200, c=COLORS['secondary'], marker='*', linewidths=0.5,
               edgecolors='k', zorder=20, label='Virginia 2011 Mw 5.8')
    ax.legend(loc='lower left', frameon=True, framealpha=0.9, fontsize=8)
    ax.set_title('Example 00: Global 1D PREM -- Station Network', fontsize=11)
    save_fig(fig, os.path.join(FIG_DIR, 'ex00_station_map'))
except ImportError:
    print('  Cartopy not available, skipping map.')

# ==== Figure 2: Seismograms at IU.ANMO =================================== #
print('Ex00 Fig 2: Seismograms')
station_key = 'IU.ANMO'
gsn_dir = os.path.join(INPUT, 'stations', 'global_seismic_network_GSN')
time = np.loadtxt(os.path.join(gsn_dir, 'data_time.ascii'))
disp = np.loadtxt(os.path.join(gsn_dir, f'{station_key}.ascii'))

# epicentral distance source -> IU.ANMO
anmo = {f'{r[1]}.{r[0]}': (float(r[2]), float(r[3])) for r in info_GSN}.get(station_key)
anmo_delta = epicentral_distance_deg(ev_lat, ev_lon, *anmo) if anmo else None

fig, axes = plt.subplots(3, 1, figsize=(7, 4.5), sharex=True, sharey=True)
for ich, ch in enumerate('RTZ'):
    ax = axes[ich]
    ax.plot(time, disp[:, ich], color=CHANNEL_COLORS[ch], lw=0.8)
    ax.set_ylabel(f'{CHANNEL_LABELS[ch]}\n(m)', fontsize=9)
    ax.tick_params(axis='both', which='both', direction='in')
    if ich < 2:
        ax.tick_params(labelbottom=False)
_dtxt = f' ($\\Delta$ = {anmo_delta:.1f}$^\\circ$)' if anmo_delta else ''
axes[0].set_title(f'Station {station_key}{_dtxt} -- 1D PREM', fontsize=11)
axes[-1].set_xlabel('Time after source origin (s)')
axes[0].set_xlim(time[0], time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'ex00_seismograms_ANMO'))

# ==== Figure 3: USArray wavefield time-lapse ============================= #
print('Ex00 Fig 3: USArray wavefield time-lapse')
info_US = np.loadtxt(os.path.join(INPUT_DIR, 'US_TA.txt'), dtype=str, skiprows=3)
nstation = len(info_US)
stlatlon = {}
for ist in range(nstation):
    key = info_US[ist, 1] + '.' + info_US[ist, 0]
    stlatlon[key] = [float(info_US[ist, 2]), float(info_US[ist, 3])]

us_dir = os.path.join(INPUT, 'stations', 'USArray_transportable')
rank_info = np.loadtxt(os.path.join(us_dir, 'rank_station.info'), dtype=str, skiprows=1)
rank_dict = {}
coords_order = []
for item in rank_info:
    rank, stkey = item[0], item[1]
    rank_dict.setdefault(rank, []).append(stkey)
    coords_order.append(stlatlon[stkey])
coords_order = np.array(coords_order)

time_us = np.loadtxt(os.path.join(us_dir, 'data_time.ascii'))
ntime = len(time_us)

# Snapshot window: span the part of the record where the wave is over the array.
# Default to (5 s, end of record); the helper picks evenly-spaced snapshots and
# auto-scales the colour limits to the 99.5th percentile of |amplitude|.
t_lo = max(5.0, float(time_us[0]))
t_hi = float(time_us[-1])
for channel, label in [('U3', 'Vertical displacement (m)'),
                       ('R3', 'Vertical rotation (rad)')]:
    data = np.ndarray((ntime, nstation))
    pos = 0
    for rank in rank_dict:
        d = np.loadtxt(os.path.join(us_dir, f'dir_rank{rank}', f'{channel}.ascii'))
        data[:, pos:pos + len(rank_dict[rank])] = d
        pos += len(rank_dict[rank])
    wavefield_timelapse_map(
        coords_order[:, 1], coords_order[:, 0], data, time_us,
        os.path.join(FIG_DIR, f'ex00_USArray_{channel}_timelapse'),
        clabel=label, title=f'1D PREM -- USArray {channel} wavefield',
        t_range=(t_lo, t_hi))

print('Ex00 figures complete.')
