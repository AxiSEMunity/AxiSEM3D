#!/usr/bin/env python3
"""
Example 01: Global 3D S362ANI — publication figures.
  Fig 1: 1D vs 3D seismogram comparison at IU.ANMO
  Fig 2: USArray U3 (vertical displacement) wavefield time-lapse (3D)
  Fig 3: source-receiver configuration map

All inputs are read RELATIVE to this script, so the example is self-contained:
run the simulation first (see README), then run this script. No absolute paths.
"""
import sys, os
import numpy as np

# ---- locate this script and the shared plotting helper ----
HERE = os.path.dirname(os.path.abspath(__file__))
# pub_style.py lives one level up, at examples/pub_style.py
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import (apply_style, save_fig, CHANNEL_COLORS, CHANNEL_LABELS,
                       epicentral_distance_deg, wavefield_timelapse_map,
                       source_receiver_map)
import matplotlib.pyplot as plt  # noqa: E402  (pub_style sets the Agg backend)

apply_style()

# ---- INPUT: where to read the simulation output from ----
# This is the ./output/ directory produced by running this example
# (mpirun -np <N> axisem3d --input input --output output). The station ascii
# files live under <INPUT>/stations/<group>/ exactly as written by the solver.
INPUT = os.path.join(HERE, 'output')
# 1D PREM run (Example 00) used only for the optional dashed comparison trace.
INPUT_1D = os.path.join(HERE, '..', '00_global_1D', 'output')
# where figures are written
FIG_DIR = os.path.join(HERE, 'figures')

# event location (Virginia 2011), used for the epicentral distance + map
EV_LAT, EV_LON = 37.91, -77.93

# ==== Figure 1: 1D vs 3D seismograms at IU.ANMO ====
print('Ex01 Fig 1: 1D vs 3D seismogram comparison')
station_key = 'IU.ANMO'
gsn_3d = os.path.join(INPUT, 'stations', 'global_seismic_network_GSN')
time = np.loadtxt(os.path.join(gsn_3d, 'data_time.ascii'))
disp3D = np.loadtxt(os.path.join(gsn_3d, f'{station_key}.ascii'))

# optional 1D overlay (only if Example 00 has been run)
gsn_1d = os.path.join(INPUT_1D, 'stations', 'global_seismic_network_GSN')
compare_1D = os.path.isfile(os.path.join(gsn_1d, f'{station_key}.ascii'))
if compare_1D:
    disp1D = np.loadtxt(os.path.join(gsn_1d, f'{station_key}.ascii'))

fig, axes = plt.subplots(3, 1, figsize=(7, 5), sharex=True)
for ich, ch in enumerate('RTZ'):
    ax = axes[ich]
    ax.plot(time, disp3D[:, ich] * 1e6, color=CHANNEL_COLORS[ch], lw=0.8,
            label='3D S362ANI')
    if compare_1D:
        ax.plot(time, disp1D[:, ich] * 1e6, color='#636363', lw=0.6,
                ls='--', label='1D PREM')
    ax.set_ylabel(f'{CHANNEL_LABELS[ch]}\n($\\mu$m)', fontsize=9)
    ax.tick_params(direction='in')
    if ich < 2:
        ax.tick_params(labelbottom=False)
axes[0].legend(loc='upper right', fontsize=8, ncol=2)
_d = epicentral_distance_deg(EV_LAT, EV_LON, 34.9462, -106.4567)  # -> IU.ANMO
axes[0].set_title(f'Station {station_key} ($\\Delta$ = {_d:.1f}$^\\circ$) '
                  f'— 1D vs 3D comparison', fontsize=11)
axes[-1].set_xlabel('Time after source origin (s)')
axes[0].set_xlim(time[0], time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'ex01_1D_vs_3D_ANMO'))

# ==== Figure 2: USArray U3 wavefield time-lapse (3D) ====
print('Ex01 Fig 2: USArray wavefield time-lapse')
info_US = np.loadtxt(os.path.join(HERE, 'input', 'US_TA.txt'), dtype=str, skiprows=3)
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

channel = 'U3'
data = np.ndarray((ntime, nstation))
pos = 0
for rank in rank_dict:
    d = np.loadtxt(os.path.join(us_dir, f'dir_rank{rank}', f'{channel}.ascii'))
    data[:, pos:pos + len(rank_dict[rank])] = d
    pos += len(rank_dict[rank])

# Snapshots evenly spaced over the propagation window, shared RdBu_r colour scale.
wavefield_timelapse_map(
    coords_order[:, 1], coords_order[:, 0], data, time_us,
    os.path.join(FIG_DIR, f'ex01_USArray_{channel}_timelapse'),
    clabel='Vertical displacement (m)',
    title=f'3D S362ANI — USArray {channel} wavefield')

# ==== Figure 3: source-receiver configuration map ====
print('Ex01 Fig 3: source-receiver map')
_recv = [(float(r[2]), float(r[3])) for r in info_US]
_sub = _recv[::max(1, len(_recv) // 120)]
source_receiver_map((EV_LAT, EV_LON), _sub,
                    os.path.join(FIG_DIR, 'ex01_source_receiver_map'),
                    title='Virginia 2011 + USArray TA — source-receiver geometry',
                    coastlines=True, recv_label=f'USArray TA (n={len(_recv)})')

print('Ex01 figures complete.')
