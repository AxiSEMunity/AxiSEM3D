#!/usr/bin/env python3
"""
Example 09: Enceladus icesheet (regional) — publication figures.

  Fig 1: Vertical-component record section, traces sorted by epicentral distance.
  Fig 2: Station map (lat/lon) with the source marked.
  Fig 3: Three-component (R/T/Z) seismogram at a representative station.
  Fig 4: Combined overview — interpolated surface wavefield snapshots + record
         section (the surface field is interpolated across the real station grid,
         requires SciPy; cartopy is used for the map projection if available).

Inputs are located RELATIVE to this script: the station ASCII produced by running
the example (output/stations/JW_stations/) plus the station-info file shipped in
input/3Dstations.txt. The shared plotting helper pub_style.py is imported from the
top-level examples/ directory (one level above this example).

Run the simulation first (see README), then:
    python make_figures.py
Figures (PNG + PDF) are written to figures/.
"""
import os
import sys
import glob

import numpy as np

import matplotlib
matplotlib.use('Agg')  # headless / no display required
import matplotlib.pyplot as plt

# locate everything relative to this file so the script is portable
HERE = os.path.dirname(os.path.abspath(__file__))
# shared style helper lives in the examples root (one level up); import relatively
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import (apply_style, save_fig, COLORS,  # noqa: E402
                       CHANNEL_COLORS, CHANNEL_LABELS)

apply_style()

# -- input/output locations (all relative to this example) --------------------
# AxiSEM3D run output produced by this example (override if you ran elsewhere):
INPUT = os.path.join(HERE, 'output')
# station group name from input/inparam.output.yaml (list_of_station_groups):
STATION_GROUP = 'JW_stations'
ST_DIR = os.path.join(INPUT, 'stations', STATION_GROUP)
STATION_FILE = os.path.join(HERE, 'input', '3Dstations.txt')
FIG_DIR = os.path.join(HERE, 'figures')

# source location (from input/inparam.source.yaml): lat/lon in degrees, depth 2 km
SRC_LAT, SRC_LON = 0.0, 0.0

# -- read station coordinates -------------------------------------------------
# 3Dstations.txt columns: NAME NETWORK latitude longitude [useless] depth
st_info = np.loadtxt(STATION_FILE, dtype=str, skiprows=2)
st_lats = st_info[:, 2].astype(float)
st_lons = st_info[:, 3].astype(float)
# map "network.name" -> (lat, lon)
name_to_latlon = {f'{row[1]}.{row[0]}': (float(row[2]), float(row[3]))
                  for row in st_info}


def epicentral_distance_deg(lat, lon):
    """Great-circle angle from the source to (lat, lon), in degrees."""
    la, lo = np.radians(lat), np.radians(lon)
    sla, slo = np.radians(SRC_LAT), np.radians(SRC_LON)
    cosd = np.sin(sla) * np.sin(la) + np.cos(sla) * np.cos(la) * np.cos(lo - slo)
    return np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0)))


# -- read the run output ------------------------------------------------------
if not os.path.isdir(ST_DIR):
    raise SystemExit(
        f'Station output not found in {ST_DIR}.\n'
        'Run the simulation first (see README "Reproduce" section), e.g.\n'
        '    mpirun -np 6 axisem3d --input input --output output')

time = np.loadtxt(os.path.join(ST_DIR, 'data_time.ascii'))

st_files = sorted(glob.glob(os.path.join(ST_DIR, '*.ascii')))
st_files = [f for f in st_files if 'data_time' not in f]
data_all = {os.path.basename(f).replace('.ascii', ''): np.loadtxt(f)
            for f in st_files}

# ==== Figure 1: Record section sorted by epicentral distance =================
print('Ex09 Fig 1: Enceladus record section')
ch_idx = 2  # vertical component (RTZ -> column 2 = Z)

# pair each station with its epicentral distance, drop any without coords
records = []
for name, d in data_all.items():
    if name in name_to_latlon:
        lat, lon = name_to_latlon[name]
        records.append((epicentral_distance_deg(lat, lon), name, d))
records.sort(key=lambda r: r[0])
dists = np.array([r[0] for r in records])

# amplitude gain (in degrees) for the normalized wiggles
dist_span = max(dists.max() - dists.min(), 1.0)
gain = dist_span / 22.0

fig, ax = plt.subplots(figsize=(7, 8))
for dist, name, d in records:
    trace = d[:, ch_idx] if d.ndim > 1 else d
    norm_val = max(np.max(np.abs(trace)), 1e-30)
    ax.plot(time, trace / norm_val * gain + dist, color=COLORS['mono'], lw=0.4)

ax.set_xlabel('Time after source origin (s)')
ax.set_ylabel('Epicentral distance ($^\\circ$)')
ax.set_xlim(time[0], time[-1])
ax.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
ax.set_title('Enceladus icesheet — Vertical component record section', fontsize=15)
save_fig(fig, os.path.join(FIG_DIR, 'ex09_enceladus_record_section'))

# ==== Figure 2: Station map with source ======================================
print('Ex09 Fig 2: Station map')
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(st_lons, st_lats, s=15, c=COLORS['primary'], marker='^',
           edgecolors='k', linewidths=0.3, zorder=10, label='Stations')
ax.scatter([SRC_LON], [SRC_LAT], s=320, c=COLORS['secondary'], marker='*',
           edgecolors='k', linewidths=0.6, zorder=20, label='Source')
ax.set_xlabel('Longitude ($^\\circ$)')
ax.set_ylabel('Latitude ($^\\circ$)')
ax.set_title('Enceladus station network', fontsize=15)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3, linewidth=0.3)
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
save_fig(fig, os.path.join(FIG_DIR, 'ex09_enceladus_station_map'))

# ==== Figure 3: Three-component seismograms ==================================
print('Ex09 Fig 3: Three-component seismograms')


def comp_onset(trace, frac=0.05):
    """First time a component exceeds frac of its OWN peak (per-component,
    so a weak-but-early transverse arrival is still detected)."""
    m = np.max(np.abs(trace))
    if m <= 0:
        return None
    return time[int(np.argmax(np.abs(trace) > frac * m))]


# pick the station where the transverse arrives earliest relative to radial
best = None
for dist, name, dd in records:
    if dd.ndim < 2 or dd.shape[1] < 3:
        continue
    oR, oT = comp_onset(dd[:, 0]), comp_onset(dd[:, 1])
    if oR is None or oT is None:
        continue
    key = oT - oR  # most negative => transverse leads radial
    if best is None or key < best[0]:
        best = (key, dist, name, dd)
rep_dist, rep_name, d = (best[1], best[2], best[3]) if best else \
    (records[len(records) // 3][0], records[len(records) // 3][1],
     data_all[records[len(records) // 3][1]])

# window must start before the onset of EVERY component (use the earliest)
onsets = [comp_onset(d[:, ich]) for ich in range(min(3, d.shape[1]))]
onsets = [o for o in onsets if o is not None]
onset_t = min(onsets) if onsets else time[0]
x0 = max(time[0], onset_t - 5.0)

fig, axes = plt.subplots(3, 1, figsize=(7, 4.5), sharex=True, sharey=True)
for ich, ch in enumerate('RTZ'):
    ax = axes[ich]
    trace = d[:, ich] if d.ndim > 1 and ich < d.shape[1] else np.zeros_like(time)
    ax.plot(time, trace, color=CHANNEL_COLORS[ch], lw=0.8)
    ax.set_ylabel(f'{CHANNEL_LABELS[ch]}\n(m)', fontsize=13)
    ax.tick_params(direction='in')
    if ich < 2:
        ax.tick_params(labelbottom=False)

axes[0].set_title(f'Enceladus seismogram — Station {rep_name} '
                  f'($\\Delta$ = {rep_dist:.1f}$^\\circ$)', fontsize=15)
axes[-1].set_xlabel('Time after source origin (s)')
axes[0].set_xlim(x0, time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'ex09_enceladus_seismograms'))

# ==== Figure 4: Combined overview — surface wavefield + record section =======
# Unlike Mars/Moon (1D AXISYMMETRIC, where the surface field is U(Delta) and is built
# by interpolating one meridian profile by epicentral distance), this is a genuinely
# 3D regional icesheet model (moment-tensor source, Nr = 5) whose receivers form a
# real 2D lat-lon GRID. We therefore do NOT collapse to U(Delta). Instead we interpolate
# the ACTUAL recorded vertical displacement across the stations onto a lat/lon grid at
# each snapshot time (scipy griddata, linear, automatically NaN/masked outside the
# station convex hull). This shows the true, non-axisymmetric surface wavefield exactly
# as sampled by the network -- no azimuthal structure is invented (stations at the same
# Delta but different azimuth keep their own amplitudes). The longitudinal sampling is
# coarse, so the field is faceted in azimuth; this is the real resolution of the array,
# not a smoothing choice. This figure needs SciPy; cartopy is used for the projection
# if available, otherwise a plain lat/lon grid is drawn.
print('Ex09 Fig 4: combined overview (interpolated surface wavefield + record section)')
try:
    from scipy.interpolate import griddata
    have_scipy = True
except ImportError:
    have_scipy = False
    print('  (skipped: SciPy not available; install scipy to render this figure)')

if have_scipy:
    import matplotlib as mpl
    import matplotlib.gridspec as gridspec

    try:
        import cartopy.crs as ccrs
        have_cartopy = True
    except ImportError:
        have_cartopy = False

    # station coordinates + vertical traces (stations with coords and a Z column)
    ov_lon, ov_lat, ov_uz = [], [], []
    for dist, name, dd in records:
        if dd.ndim > 1 and dd.shape[1] > ch_idx:
            la, lo = name_to_latlon[name]
            ov_lon.append(lo)
            ov_lat.append(la)
            ov_uz.append(dd[:, ch_idx])
    ov_lon = np.array(ov_lon)
    ov_lat = np.array(ov_lat)
    ov_uz = np.array(ov_uz)                                  # (nstation, ntime)
    pts = np.column_stack([ov_lon, ov_lat])

    # The vertical amplitude decays by orders of magnitude from the near-source station
    # to Delta ~ 60 deg, so an absolute-amplitude snapshot is saturated by the source
    # blob and the travelling wavefront is invisible. We therefore normalise EACH station
    # to its own peak before interpolating -- identical to the per-trace normalisation in
    # the record section -- so the panels show the PHASE/ARRIVAL of the wavefront across
    # the network on a common +/-1 scale (a display choice, not invented data).
    peak = np.maximum(np.max(np.abs(ov_uz), axis=1, keepdims=True), 1e-30)
    ov_uz_n = ov_uz / peak                                   # (nstation, ntime) in [-1, 1]

    # regular lat/lon grid spanning the station patch (small margin beyond the array)
    lon_lo, lon_hi = ov_lon.min() - 3, ov_lon.max() + 3
    lat_lo, lat_hi = ov_lat.min() - 3, ov_lat.max() + 3
    glon = np.linspace(lon_lo, lon_hi, 240)
    glat = np.linspace(lat_lo, lat_hi, 240)
    GLON, GLAT = np.meshgrid(glon, glat)

    # snapshot times: span from the first arrival at the nearest station to the first
    # arrival at the farthest, so the panels track the wavefront sweeping out to the edge.
    def _onset(tr, frac=0.05):
        m = np.max(np.abs(tr))
        return time[int(np.argmax(np.abs(tr) > frac * m))] if m > 0 else time[-1]

    arr = np.array([_onset(tr) for tr in ov_uz])
    t0 = float(arr.min())
    t1 = float(min(arr.max() + 0.10 * (arr.max() - arr.min()), time[-1]))
    six_t = np.linspace(t0, t1, 6)
    SR, SC = 3, 2

    clon = float(0.5 * (lon_lo + lon_hi))
    clat = float(0.5 * (lat_lo + lat_hi))
    snap_kw = ({'projection': ccrs.Orthographic(clon, clat)} if have_cartopy else {})

    fig = plt.figure(figsize=(11, 8.5))
    fig.set_layout_engine('none')
    outer = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1.0], wspace=0.16, figure=fig)
    gs_snap = gridspec.GridSpecFromSubplotSpec(SR, SC, subplot_spec=outer[0],
                                               wspace=0.10, hspace=0.24)

    for k, t in enumerate(six_t):
        ax = fig.add_subplot(gs_snap[k // SC, k % SC], **snap_kw)
        it = int(np.argmin(np.abs(time - t)))
        field = griddata(pts, ov_uz_n[:, it], (GLON, GLAT), method='linear')  # NaN outside hull
        field = np.ma.masked_invalid(field)
        vmax = 0.9                                            # fixed: traces are normalised
        if have_cartopy:
            ax.set_extent([lon_lo, lon_hi, lat_lo, lat_hi], crs=ccrs.PlateCarree())
            ax.pcolormesh(glon, glat, field, transform=ccrs.PlateCarree(),
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax, rasterized=True,
                          shading='auto')
            ax.scatter(ov_lon, ov_lat, s=7, c='k', transform=ccrs.PlateCarree(),
                       zorder=8, linewidths=0)
            ax.plot(SRC_LON, SRC_LAT, '*', ms=16, color='k',
                    transform=ccrs.PlateCarree(), zorder=10)
            ax.plot(SRC_LON, SRC_LAT, '*', ms=10, color=COLORS['secondary'],
                    transform=ccrs.PlateCarree(), zorder=11)
            # lat/lon graticule so the spatial scale is explicit; label only the outer
            # edges (left column -> lat, bottom row -> lon) to keep the small panels clean
            gl = ax.gridlines(draw_labels=True, linewidth=0.6, color='0.55', alpha=0.6,
                              xlocs=np.arange(-30, 61, 15), ylocs=np.arange(-60, 61, 30))
            gl.top_labels = gl.right_labels = False
            gl.bottom_labels = (k // SC == SR - 1)
            gl.left_labels = (k % SC == 0)
            gl.xlabel_style = {'size': 13}
            gl.ylabel_style = {'size': 13}
        else:
            ax.pcolormesh(glon, glat, field, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                          rasterized=True, shading='auto')
            ax.scatter(ov_lon, ov_lat, s=7, c='k', zorder=8)
            ax.plot(SRC_LON, SRC_LAT, '*', ms=14, color=COLORS['secondary'], zorder=10)
            ax.set_aspect('equal')
            ax.set_xticks(np.arange(-30, 61, 15))
            ax.set_yticks(np.arange(-60, 61, 30))
            ax.tick_params(labelsize=13)
            ax.grid(True, lw=0.6, color='0.6', alpha=0.5)
            if k % SC == 0:
                ax.set_ylabel('lat ($^\\circ$)', fontsize=15)
            if k // SC == SR - 1:
                ax.set_xlabel('lon ($^\\circ$)', fontsize=15)
        ax.set_title(f't = {time[it]:.0f} s', fontsize=17)

    # right: vertical-component record section (same gain/style as Fig 1)
    axR = fig.add_subplot(outer[1])
    for dist, name, dd in records:
        trace = dd[:, ch_idx] if dd.ndim > 1 else dd
        norm_val = max(np.max(np.abs(trace)), 1e-30)
        axR.plot(time, trace / norm_val * gain + dist, color=COLORS['mono'], lw=0.4)
    axR.set_xlabel('Time after source origin (s)', fontsize=17)
    axR.set_ylabel('Epicentral distance ($^\\circ$)', fontsize=17)
    axR.set_xlim(time[0], time[-1])
    axR.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
    axR.tick_params(labelsize=15)
    axR.set_title('Vertical-component record section', fontsize=17)

    fig.suptitle('Enceladus icesheet — interpolated surface wavefield and record section',
                 fontsize=20)
    fig.subplots_adjust(left=0.04, right=0.97, top=0.90, bottom=0.13)
    # wavefield amplitude scale (the surface field is normalised per station to its peak)
    cax = fig.add_axes([0.07, 0.045, 0.33, 0.014])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                      cax=cax, orientation='horizontal')
    cb.set_label('$U_z$ / station peak  (each station normalised)', fontsize=13)
    cb.set_ticks([-1, 0, 1])
    cb.ax.tick_params(labelsize=8)
    save_fig(fig, os.path.join(FIG_DIR, 'ex09_enceladus_overview'), tight=False)

print('Ex09 figures complete.')
