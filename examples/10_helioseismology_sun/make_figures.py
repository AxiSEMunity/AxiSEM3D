#!/usr/bin/env python3
"""
Helioseismology (Sun) AxiSEM3D simulation — publication figures.
  Fig 1: receiver map (polar view around the pole source)
  Fig 2: three-component displacement seismograms at one receiver
  Fig 3: vertical-displacement record section vs. epicentral distance

Reads ASCII_STATION output (one file per receiver, columns = R,T,Z) that this
example produces in output/stations/stationsforplotting/, plus the source and
receiver lists from input/. Imports only the shared pub_style helper.

Run from anywhere with a Python that has numpy + pyyaml + matplotlib:
    python make_figures.py
"""
import sys, os
import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))

# pub_style.py is the shared plotting-style helper staged in the examples root
# (one level above this example directory). Add that directory to sys.path
# relative to this file so the import works wherever the repo is checked out.
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import apply_style, save_fig, COLORS, CHANNEL_COLORS, CHANNEL_LABELS
import matplotlib.pyplot as plt

apply_style()

# INPUT: directory holding this example's AxiSEM3D run output (its ./output/).
INPUT = os.path.join(HERE, 'output')
STA_DIR = os.path.join(INPUT, 'stations', 'stationsforplotting')
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)


def gc_distance_deg(lat1, lon1, lat2, lon2):
    """Great-circle distance in degrees."""
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dl = np.radians(lon2 - lon1)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


# ---- source location ----
with open(os.path.join(HERE, 'input', 'inparam.source.yaml')) as f:
    src = yaml.safe_load(f)
s0 = list(src['list_of_sources'][0].values())[0]
ev_lat, ev_lon = s0['location']['latitude_longitude']

# ---- read receiver list ----
rows = []
with open(os.path.join(HERE, 'input', 'stations_for_plotting.txt')) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        c = line.split()
        rows.append((c[0], c[1], float(c[2]), float(c[3])))  # name, net, lat, lon
stations = {f'{net}.{nam}': (lat, lon) for nam, net, lat, lon in rows}

# ---- load traces ----
time = np.loadtxt(os.path.join(STA_DIR, 'data_time.ascii'))
traces = {}
for key, (lat, lon) in stations.items():
    fpath = os.path.join(STA_DIR, f'{key}.ascii')
    if not os.path.exists(fpath):
        continue
    d = np.loadtxt(fpath)
    if d.ndim == 1:
        d = d[:, None]
    traces[key] = (lat, lon, d)
print(f'Loaded {len(traces)} / {len(stations)} receiver files; {len(time)} time samples')

# ==== Figure 1: receiver map (polar, source-centred) ====
print('Sun Fig 1: receiver map')
fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={'projection': 'polar'})
ax.set_theta_zero_location('N')
for key, (lat, lon, _) in traces.items():
    dist = gc_distance_deg(ev_lat, ev_lon, lat, lon)
    ax.plot(np.radians(lon), dist, '^', ms=5, color=COLORS['primary'],
            markeredgecolor='k', markeredgewidth=0.3)
ax.plot(np.radians(ev_lon), 0, '*', ms=18, color=COLORS['secondary'],
        markeredgecolor='k', markeredgewidth=0.5, zorder=20)
ax.set_rlabel_position(135)
ax.set_title('Sun: pressure source (star) and receivers\n(radius = angular distance, deg)',
             fontsize=10, pad=18)
save_fig(fig, os.path.join(FIG_DIR, 'sun_receiver_map'))

# ==== Figure 2: three-component seismograms at one receiver ====
print('Sun Fig 2: seismograms')
# pick a mid-distance receiver
keys_by_dist = sorted(traces, key=lambda k: gc_distance_deg(ev_lat, ev_lon, traces[k][0], traces[k][1]))
pick = keys_by_dist[len(keys_by_dist) // 2]
lat, lon, d = traces[pick]
dpick = gc_distance_deg(ev_lat, ev_lon, lat, lon)
ncomp = min(3, d.shape[1])
fig, axes = plt.subplots(ncomp, 1, figsize=(7, 4.5), sharex=True)
if ncomp == 1:
    axes = [axes]
for ich in range(ncomp):
    ch = 'RTZ'[ich]
    axes[ich].plot(time, d[:, ich], color=CHANNEL_COLORS[ch], lw=0.7)
    axes[ich].set_ylabel(f'{CHANNEL_LABELS[ch]}\n$U$ (m)', fontsize=9)
    axes[ich].tick_params(direction='in')
axes[0].set_title(f'Sun receiver {pick}  (Δ = {dpick:.0f}°)', fontsize=11)
axes[-1].set_xlabel('Time after source origin (s)')
axes[0].set_xlim(time[0], time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'sun_seismograms'))

# ==== Figure 3: vertical-displacement record section ====
print('Sun Fig 3: record section')
# take one azimuth column if present (network 'a0'); else all receivers
sel = [(gc_distance_deg(ev_lat, ev_lon, lat, lon), key, d)
       for key, (lat, lon, d) in traces.items() if key.startswith('a0.')]
if not sel:
    sel = [(gc_distance_deg(ev_lat, ev_lon, lat, lon), key, d)
           for key, (lat, lon, d) in traces.items()]
sel.sort()
fig, ax = plt.subplots(figsize=(7, 6))
zc = min(2, sel[0][2].shape[1] - 1)  # vertical channel index
# normalise each trace, offset by distance
for dist, key, d in sel:
    tr = d[:, zc].astype(float)
    norm = np.max(np.abs(tr))
    if norm > 0:
        tr = tr / norm * 6.0  # amplitude in degrees of offset
    ax.plot(time, dist + tr, color=COLORS['mono'], lw=0.5)
ax.set_xlabel('Time after source origin (s)')
ax.set_ylabel('Epicentral distance (deg)')
ax.set_xlim(time[0], time[-1])
ax.invert_yaxis()
ax.set_title('Sun: vertical-displacement record section', fontsize=11)
save_fig(fig, os.path.join(FIG_DIR, 'sun_record_section'))

print('Sun figures complete.')
