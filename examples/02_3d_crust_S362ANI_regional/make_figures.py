#!/usr/bin/env python3
"""
Example 02: Regional 3-D crust (CRUST1.0) on S362ANI — publication figures.

Reproduces, directly from this example's own simulation output:
  Fig 1: 1D vs 3D crust record section (vertical displacement, by distance)
  Fig 2: USArray U3 wavefield time-lapse, 1D crust
  Fig 3: USArray U3 wavefield time-lapse, 3D crust (same colour scale as Fig 2)
  Fig 4: source-receiver geometry map (Virginia 2011 + USArray TA)

Inputs are the two runs assembled and launched per the README (the
simu_with_1d_crust / simu_with_3d_crust cases): each writes ASCII_CHANNEL
station output under
    simu_with_<case>_crust/output/stations/USArray_transportable/
(per-rank dir_rank*/U3.ascii, rank_station.info, data_time.ascii). Station
coordinates are read from the shipped input_share/US_TA.txt.

Run after both simulations have finished:
    python make_figures.py
Figures are written (PDF + PNG) into ./figures/.
"""
import sys
import os
import numpy as np

# Locate inputs relative to this script so the example is portable.
HERE = os.path.dirname(os.path.abspath(__file__))
# The shared pub_style helper lives one level up, at examples/pub_style.py.
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import (apply_style, save_fig, COLORS,
                       epicentral_distance_deg, wavefield_timelapse_map,
                       source_receiver_map)
import matplotlib.pyplot as plt  # noqa: E402  (pub_style sets the Agg backend)

apply_style()

# ---- INPUT: the two simulation outputs this example produces -----------------
# Each run writes its station group "USArray_transportable" here. Override these
# if you ran the solver with a different --output directory.
STATION_GROUP = 'USArray_transportable'
INPUT_1D = os.path.join(HERE, 'simu_with_1d_crust', 'output', 'stations', STATION_GROUP)
INPUT_3D = os.path.join(HERE, 'simu_with_3d_crust', 'output', 'stations', STATION_GROUP)

# Receiver coordinates (shipped with the example).
US_TA_FILE = os.path.join(HERE, 'input_share', 'US_TA.txt')

# Where figures are written.
FIG_DIR = os.path.join(HERE, 'figures')

# 2011 Virginia event (matches input_share/inparam.source.yaml).
SRC_LAT, SRC_LON = 37.91, -77.93

# ---- station table -----------------------------------------------------------
# US_TA.txt columns: NAME NETWORK LATITUDE LONGITUDE [useless] [depth]
info_US = np.loadtxt(US_TA_FILE, dtype=str, skiprows=3)
nstation = len(info_US)
stlatlon = {}
for ist in range(nstation):
    key = info_US[ist, 1] + '.' + info_US[ist, 0]   # NETWORK.NAME
    stlatlon[key] = [float(info_US[ist, 2]), float(info_US[ist, 3])]


def read_us_ta(us_ta_dir, channel='U3'):
    """Read an ASCII_CHANNEL station group into (time, data, coords, names).

    data has shape (ntime, nstation); coords is (nstation, 2) = (lat, lon)."""
    rank_info = np.loadtxt(os.path.join(us_ta_dir, 'rank_station.info'),
                           dtype=str, skiprows=1)
    rank_dict = {}
    coords_order = []
    names_order = []
    for item in rank_info:
        rank, stkey = item[0], item[1]
        rank_dict.setdefault(rank, []).append(stkey)
        coords_order.append(stlatlon[stkey])
        names_order.append(stkey)
    coords_order = np.array(coords_order)
    time = np.loadtxt(os.path.join(us_ta_dir, 'data_time.ascii'))
    ntime = len(time)
    data = np.ndarray((ntime, len(names_order)))
    pos = 0
    for rank in rank_dict:
        d = np.loadtxt(os.path.join(us_ta_dir, f'dir_rank{rank}', f'{channel}.ascii'))
        data[:, pos:pos + len(rank_dict[rank])] = d
        pos += len(rank_dict[rank])
    return time, data, coords_order, names_order


if not os.path.isdir(INPUT_1D) or not os.path.isdir(INPUT_3D):
    print('Ex02: simulation output not found.')
    print(f'  expected 1D crust output: {INPUT_1D}')
    print(f'  expected 3D crust output: {INPUT_3D}')
    print('  Run the 1D-crust and 3D-crust cases first (see README, "Run").')
    sys.exit(0)

time_1d, data_1d, coords_1d, names_1d = read_us_ta(INPUT_1D)
time_3d, data_3d, coords_3d, names_3d = read_us_ta(INPUT_3D)

# ==== Figure 1: Record section comparison =====================================
print('Ex02 Fig 1: 1D vs 3D crust record section')
# epicentral distance of every station, then select a subset spread evenly
# across the distance range so the record section shows clean moveout
dist_all = np.array([epicentral_distance_deg(SRC_LAT, SRC_LON, lat, lon)
                     for lat, lon in coords_1d])
order = np.argsort(dist_all)
n_seis = min(60, len(order))
sel = order[np.linspace(0, len(order) - 1, n_seis).astype(int)]
dsel = dist_all[sel]
span = max(dsel.max() - dsel.min(), 1.0)
gain = span / 28.0   # amplitude of the normalized wiggles, in degrees

fig, ax = plt.subplots(figsize=(7, 9))
for j, (ind, dist) in enumerate(zip(sel, dsel)):
    norm_val = max(np.max(np.abs(data_1d[:, ind])), 1e-30)
    ax.plot(time_1d, data_1d[:, ind] / norm_val * gain + dist,
            color=COLORS['primary'], lw=0.5, label='1D crust' if j == 0 else None)
    ax.plot(time_3d, data_3d[:, ind] / norm_val * gain + dist,
            color=COLORS['secondary'], lw=0.4, label='3D crust' if j == 0 else None)

ax.set_xlabel('Time after source origin (s)')
ax.set_ylabel('Epicentral distance ($^\\circ$)')
ax.set_xlim(time_1d[0], time_1d[-1])
ax.set_ylim(dsel.min() - gain * 1.5, dsel.max() + gain * 1.5)
ax.legend(loc='upper right', fontsize=8)
ax.set_title('1D vs 3D crustal model — Vertical displacement', fontsize=11)
save_fig(fig, os.path.join(FIG_DIR, 'ex02_1D_vs_3D_crust_record_section'))

# ==== Figures 2 & 3: Wavefield time-lapse (1D and 3D crust) ====================
print('Ex02 Figs 2 & 3: USArray U3 wavefield time-lapse')
# Both time-lapse grids share ONE colour scale so the 1D and 3D crust runs are
# directly comparable. Use the larger of the two runs' 99.5th-percentile |U3|
# over the chosen snapshot window.
_t_range = (5, min(1200, float(time_1d[-1])))


def _window_vmax(data, times, t_range, pct=99.5):
    target = np.linspace(t_range[0], t_range[1], 8)
    snap = np.array([int(np.argmin(np.abs(times - tt))) for tt in target])
    return float(np.percentile(np.abs(data[snap]), pct)) or 1e-30


_vmax = max(_window_vmax(data_1d, time_1d, _t_range),
            _window_vmax(data_3d, time_3d, _t_range))

for data, coords, tag, label in [(data_1d, coords_1d, '1D', '1D crust'),
                                 (data_3d, coords_3d, '3D', '3D crust')]:
    wavefield_timelapse_map(
        coords[:, 1], coords[:, 0], data, time_1d,
        os.path.join(FIG_DIR, f'ex02_USArray_{tag}_crust_timelapse'),
        clabel='Vertical displacement (m)', vmax=_vmax, t_range=_t_range,
        title=f'{label} — USArray U3 wavefield')

# ==== Figure 4: source-receiver configuration map =============================
print('Ex02 Fig 4: source-receiver map')
_recv = [(float(r[2]), float(r[3])) for r in info_US]
# thin the dense TA to ~120 markers so the map stays legible
_sub = _recv[::max(1, len(_recv) // 120)]
source_receiver_map((SRC_LAT, SRC_LON), _sub,
                    os.path.join(FIG_DIR, 'ex02_source_receiver_map'),
                    title='Virginia 2011 + USArray TA — source-receiver geometry',
                    coastlines=True, recv_label=f'USArray TA (n={len(_recv)})')

print('Ex02 figures complete.')
