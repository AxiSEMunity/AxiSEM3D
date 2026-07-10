#!/usr/bin/env python3
"""
Example 08: Mars Atmosphere (global, acoustic-seismic coupling) -- figures.

  Fig 1: ex08_mars_record_section -- vertical-component (Z) record section,
         traces normalised and offset by epicentral distance (0 deg .. 9 deg).
  Fig 2: ex08_mars_seismograms    -- three-component (R/T/Z) seismograms at a
         representative low-distance station (R/T/Z share one y-scale).
  Fig 3: ex08_source_receiver_map -- source-receiver geometry on the Mars globe
         (no Earth coastlines). Needs cartopy; skipped with a note if missing.
  Fig 4: ex08_mars_overview       -- orthographic surface-wavefield snapshots
         next to the record section. Falls back to a flat lat/lon view if
         cartopy is missing.

Note: this is a 1D axisymmetric model with the source and all stations on one
meridian, so the TRANSVERSE component is identically zero by symmetry -- its
flat panel in Fig 2 is physics, not a bug.

Inputs are located RELATIVE to this script: it reads the AxiSEM3D station
output that this example produces under
  axisem3d_mars_atm/output/stations/MARS/   (set by INPUT_DIR below)
and the station-info file axisem3d_mars_atm/input/mars.txt. Run the solver
first (see README "Reproduce"); figures are written next to this script into
figures/ as PDF + PNG.
"""
import os
import sys
import glob
import numpy as np

# this example's directory, and the examples-root (one level up) so the shared
# pub_style helper at examples/pub_style.py imports without an absolute path.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from pub_style import (apply_style, save_fig, COLORS, CHANNEL_COLORS,
                       CHANNEL_LABELS, epicentral_distance_deg,
                       source_receiver_map)
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

apply_style()

# --- inputs (all relative to this script) ---------------------------------
SIM_DIR = os.path.join(HERE, 'axisem3d_mars_atm')
# station seismograms produced by this example's run (ASCII_STATION format)
INPUT_DIR = os.path.join(SIM_DIR, 'output', 'stations', 'MARS')
# station-info file (name, network, lat, lon) used for epicentral distances
STATION_FILE = os.path.join(SIM_DIR, 'input', 'mars.txt')
FIG_DIR = os.path.join(HERE, 'figures')

# source location (inparam.source.yaml -> latitude_longitude: [0, 0])
SRC_LAT, SRC_LON = 0.0, 0.0

if not os.path.isdir(INPUT_DIR):
    sys.exit(
        f"Station output not found at:\n  {INPUT_DIR}\n"
        "Run the simulation first (see README 'Reproduce'):\n"
        "  cd axisem3d_mars_atm && mpirun -np 4 axisem3d --input input "
        "--output output")

# station info: "network.name" -> (lat, lon)
st_info = np.loadtxt(STATION_FILE, dtype=str, skiprows=3)
name2latlon = {f'{row[1]}.{row[0]}': (float(row[2]), float(row[3]))
               for row in st_info}

time = np.loadtxt(os.path.join(INPUT_DIR, 'data_time.ascii'))

st_files = sorted(glob.glob(os.path.join(INPUT_DIR, '*.ascii')))
st_files = [f for f in st_files if 'data_time' not in f]
data_all = {os.path.basename(f).replace('.ascii', ''): np.loadtxt(f)
            for f in st_files}


def dist_of(name):
    if name in name2latlon:
        return epicentral_distance_deg(SRC_LAT, SRC_LON, *name2latlon[name])
    return None


# ==== Figure 1: Record section vs epicentral distance ====
print('Ex08 Fig 1: Mars record section')
records = []
for name, d in data_all.items():
    dd = dist_of(name)
    if dd is not None:
        records.append((dd, name, d))
records.sort(key=lambda r: r[0])
dists = np.array([r[0] for r in records])

ch_idx = 2  # vertical
span = max(dists.max() - dists.min(), 1.0)
gain = span / 18.0

fig, ax = plt.subplots(figsize=(7, 5))
for dd, name, d in records:
    trace = d[:, ch_idx] if d.ndim > 1 else d
    norm_val = max(np.max(np.abs(trace)), 1e-30)
    ax.plot(time, trace / norm_val * gain + dd, color=COLORS['mono'], lw=0.7)

ax.set_xlabel('Time after source origin (s)')
ax.set_ylabel('Epicentral distance ($^\\circ$)')
ax.set_xlim(time[0], time[-1])
ax.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
ax.set_title('Mars atmosphere -- Vertical component record section', fontsize=15)
save_fig(fig, os.path.join(FIG_DIR, 'ex08_mars_record_section'))

# ==== Figure 2: Three-component seismograms ====
print('Ex08 Fig 2: Three-component seismograms')
rep_dist, rep_name, d = records[len(records) // 4]  # lower-distance station

fig, axes = plt.subplots(3, 1, figsize=(7, 4.5), sharex=True, sharey=True)
for ich, ch in enumerate('RTZ'):
    ax = axes[ich]
    trace = d[:, ich] if d.ndim > 1 and ich < d.shape[1] else np.zeros_like(time)
    ax.plot(time, trace, color=CHANNEL_COLORS[ch], lw=0.8)
    ax.set_ylabel(f'{CHANNEL_LABELS[ch]}\n(m)', fontsize=13)
    ax.tick_params(direction='in')
    if ich < 2:
        ax.tick_params(labelbottom=False)

axes[0].set_title(f'Mars seismogram -- Station {rep_name} '
                  f'($\\Delta$ = {rep_dist:.1f}$^\\circ$)', fontsize=15)
axes[-1].set_xlabel('Time after source origin (s)')
axes[0].set_xlim(time[0], time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'ex08_mars_seismograms'))

# ==== Figure 3: source-receiver configuration map (Mars globe, no coastlines) ==
print('Ex08 Fig 3: source-receiver map')
_recv = [(float(r[2]), float(r[3])) for r in st_info]
source_receiver_map((SRC_LAT, SRC_LON), _recv,
    os.path.join(FIG_DIR, 'ex08_source_receiver_map'),
    title='Mars -- source-receiver geometry', coastlines=False,
    recv_label=f'Mars stations (n={len(_recv)})')

# ==== Figure 4: combined overview (surface snapshots + record section) ========
# Mars here is a 1D AXISYMMETRIC model with an on-axis (Nr = 1) pressure source
# and all receivers on one meridian (lon 0), so the surface field is a function
# of epicentral distance Delta alone: U(Delta, t). Each map point (lat, lon)
# takes the value at Delta = arccos(cos lat * cos lon) from the source. The
# receiver line spans Delta = 0..9 deg, so the field is resolved on a near-axis
# cap; we mask everything beyond the farthest receiver (do NOT extrapolate into
# the unsampled far hemisphere) and use a zoomed orthographic so the resolved
# cap fills the frame.
#
# NOTE on the on-axis station: the pressure source sits on the mesh axis, so the
# Delta = 0 receiver (A.0) carries the near-source singularity and is ~1000x
# larger than the propagating-wave amplitudes at A.1..A.9. Including A.0 in the
# colour normalisation would saturate one central pixel and wash the travelling
# wavefront out to a flat disk. We therefore build the surface field from the
# OFF-AXIS line (A.1..A.9), the outward-propagating wave, and mask the inner
# Delta < dmin cap. A.0 is still shown in full in the record section on the right.
print('Ex08 Fig 4: combined overview (surface snapshots + record section)')

# vertical surface displacement on the receiver line, ordered by distance
recs_by_d = sorted(records, key=lambda r: r[0])
all_dist = np.array([r[0] for r in recs_by_d])               # epicentral dist (deg)
all_uz = np.array([r[2][:, ch_idx] if r[2].ndim > 1 else r[2]
                   for r in recs_by_d])                      # (nstation, ntime)
# off-axis propagating wave only (drop the Delta=0 near-source singularity)
off = all_dist > 0.5
line_dist = all_dist[off]
line_uz = all_uz[off]
dmin = float(line_dist.min())                                # inner edge of ring
dmax = float(line_dist.max())                                # resolved cap radius

# lat/lon grid -> epicentral distance from on-axis source at (0, 0)
nlat, nlon = 361, 721
lats = np.linspace(-90, 90, nlat)
lons = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lons, lats)
DELTA = np.degrees(np.arccos(np.clip(
    np.cos(np.radians(LAT)) * np.cos(np.radians(LON)), -1, 1)))
RESOLVED = (DELTA >= dmin) & (DELTA <= dmax)                 # only the sampled ring

try:
    import cartopy.crs as ccrs
    have_cartopy = True
except ImportError:
    have_cartopy = False

# snapshot times across the wave-propagation window (first arrivals -> coda)
amp_t = np.max(np.abs(line_uz), axis=0)
act = np.where(amp_t > 0.04 * max(amp_t.max(), 1e-30))[0]
t0 = time[act[0]] if act.size else time[0]
t1 = time[act[-1]] if act.size else time[-1]
six_t = np.linspace(t0, t1, 6)
SR, SC = 3, 2

# We hand-place every axis with GridSpec + add_axes, so give the figure its own
# explicit subplot margins (GridSpec bounds below) instead of letting pub_style's
# constrained_layout reflow them. Bottom leaves room for the colourbar caption.
fig = plt.figure(figsize=(11, 8.5), layout='none')
outer = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.16, figure=fig,
                          left=0.05, right=0.97, top=0.90, bottom=0.16)
gs_snap = gridspec.GridSpecFromSubplotSpec(SR, SC, subplot_spec=outer[0],
                                           wspace=0.10, hspace=0.24)
# zoomed orthographic centred on the source pole, cropped to the resolved cap.
# (set_extent in PlateCarree degrees is the reliable way to zoom a GeoAxes onto
# a small near-axis cap; set_global + set_xlim does not crop a tiny cap right.)
snap_kw = ({'projection': ccrs.Orthographic(0, 0)} if have_cartopy else {})
pad = dmax * 1.30                                            # cap half-extent (deg)

for k, t in enumerate(six_t):
    ax = fig.add_subplot(gs_snap[k // SC, k % SC], **snap_kw)
    it = int(np.argmin(np.abs(time - t)))
    field = np.interp(DELTA, line_dist, line_uz[:, it],
                      left=line_uz[0, it], right=line_uz[-1, it])
    field = np.ma.array(field, mask=~RESOLVED)              # blank unsampled area
    vmax = max(np.percentile(np.abs(field.compressed()), 99.0), 1e-30) \
        if field.count() else 1e-30
    if have_cartopy:
        ax.set_extent([-pad, pad, -pad, pad], crs=ccrs.PlateCarree())
        ax.pcolormesh(lons, lats, field, transform=ccrs.PlateCarree(),
                      cmap='RdBu_r', vmin=-vmax, vmax=vmax, rasterized=True,
                      shading='auto')
        ax.plot(0, 0, '*', ms=15, color='k', transform=ccrs.PlateCarree(),
                zorder=10)
        # epicentral-distance rings (labelled in degrees) tie the cap to the
        # record-section axis; this is the scale, so we don't add lat/lon tick
        # labels (cartopy's gridliner places them off-canvas for tiny caps).
        cs = ax.contour(LON, LAT, DELTA, levels=[3, 6, 9], colors='0.4',
                        linewidths=0.6, transform=ccrs.PlateCarree(), zorder=6)
        if k == 0:
            ax.clabel(cs, fmt=lambda v: f'{v:.0f}$^\\circ$', fontsize=11,
                      inline=True)
        ax.gridlines(draw_labels=False, linewidth=0.6, color='0.6', alpha=0.5,
                     xlocs=np.arange(-9, 10, 3), ylocs=np.arange(-9, 10, 3))
    else:
        ax.pcolormesh(lons, lats, field, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      rasterized=True, shading='auto')
        ax.contour(LON, LAT, DELTA, levels=[3, 6, 9], colors='0.4',
                   linewidths=0.6, zorder=6)
        ax.set_xlim(-pad, pad)
        ax.set_ylim(-pad, pad)
        ax.set_aspect('equal')
        ax.set_xticks(np.arange(-9, 10, 3))
        ax.set_yticks(np.arange(-9, 10, 3))
        ax.tick_params(labelsize=13)
        if k % SC == 0:
            ax.set_ylabel('lat ($^\\circ$)', fontsize=15)
        if k // SC == SR - 1:
            ax.set_xlabel('lon ($^\\circ$)', fontsize=15)
    ax.set_title(f't = {time[it]:.0f} s', fontsize=17)

# right: vertical-component record section (same gain/style as Fig 1)
axR = fig.add_subplot(outer[1])
for dd, name, d in records:
    tr = d[:, ch_idx] if d.ndim > 1 else d
    nv = max(np.max(np.abs(tr)), 1e-30)
    axR.plot(time, tr / nv * gain + dd, color=COLORS['mono'], lw=0.7)
axR.set_xlabel('Time after source origin (s)', fontsize=17)
axR.set_ylabel('Epicentral distance ($^\\circ$)', fontsize=17)
axR.set_xlim(time[0], time[-1])
axR.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
axR.tick_params(labelsize=15)
axR.set_title('Vertical-component record section', fontsize=17)

fig.suptitle(f'Mars atmosphere -- surface wavefield ({dmin:.0f}$-${dmax:.0f}'
             '$^\\circ$) and record section', fontsize=20)
# wavefield amplitude scale (each snapshot self-scaled to its own peak),
# placed beneath the snapshot column (figure coords)
cax = fig.add_axes([0.08, 0.055, 0.33, 0.014])
cb = fig.colorbar(mpl.cm.ScalarMappable(
    norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
    cax=cax, orientation='horizontal')
cb.set_label('$U_z$ / max$|U_z|$  (each snapshot normalised)', fontsize=13)
cb.set_ticks([-1, 0, 1])
cb.ax.tick_params(labelsize=8)
save_fig(fig, os.path.join(FIG_DIR, 'ex08_mars_overview'), tight=False)

print('Ex08 figures complete.')
