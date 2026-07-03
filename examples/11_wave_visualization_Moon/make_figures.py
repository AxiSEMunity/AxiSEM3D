#!/usr/bin/env python3
"""
Moon meteoroid-impact wave visualization — self-contained publication figures.
  Fig 1: lunar surface vertical-displacement snapshots (orthographic), several t
  Fig 2: surface-wave record section (epicentral distance vs. time)
  Fig 3: three-component seismogram at one lunar receiver

A vertical/isotropic on-axis impact radiates an AXISYMMETRIC wavefield, so the
surface field U(distance, t) along one meridian line of receivers fully
determines the field everywhere on the Moon: each map point (lat, lon) takes the
value at epicentral distance Delta = arccos(cos lat * cos lon) from the impact.
Reads only the compact ./data/moon_surface.npz and the local pub_style.py.
Run: <python with numpy/matplotlib(+cartopy)> make_figures.py
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style, save_fig, COLORS, CHANNEL_COLORS, CHANNEL_LABELS
import matplotlib.pyplot as plt

apply_style()
FIG_DIR = HERE
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist = d['dist']; time = d['time']; uz = d['uz']          # uz (station, time)
seis_time = d['seis_time']; seis = d['seis']; seis_dist = float(d['seis_dist'])
RES = str(d['res'])
print(f'{len(dist)} receivers, Delta = {dist.min():.1f}..{dist.max():.1f} deg, {len(time)} frames ({RES} run)')

# ==== Figure 1: lunar surface wavefield snapshots ====
print('Moon Fig 1: surface snapshots')
# 8 snapshots at finer time steps (2 x 4 grid). The orthographic view is centred
# on the impact, so it shows only the NEAR hemisphere (0-90 deg); the wavefront
# reaches the limb by ~750 s, so sample the near-side propagation window densely.
NR, NC = 2, 4
snap_t = [t for t in np.linspace(80, 740, NR * NC) if t <= time[-1]]
nlat, nlon = 180, 360
lats = np.linspace(-90, 90, nlat); lons = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lons, lats)
DELTA = np.degrees(np.arccos(np.clip(np.cos(np.radians(LAT)) * np.cos(np.radians(LON)), -1, 1)))
try:
    import cartopy.crs as ccrs
    have_cartopy = True
except ImportError:
    have_cartopy = False

n = len(snap_t)
if have_cartopy:
    fig, axes = plt.subplots(NR, NC, figsize=(2.6 * NC, 2.7 * NR),
                             subplot_kw={'projection': ccrs.Orthographic(0, 0)})
else:
    fig, axes = plt.subplots(NR, NC, figsize=(2.6 * NC, 2.4 * NR))
axes = np.atleast_1d(axes).ravel()
for ax, t in zip(axes, snap_t):
    it = int(np.argmin(np.abs(time - t)))
    field = np.interp(DELTA, dist, uz[:, it], left=uz[0, it], right=uz[-1, it])
    # per-snapshot colour scale: the wavefront amplitude decays strongly as it
    # spreads, so a single global vmax makes the late wavefronts invisible.
    # Normalise each panel to its own 99th-percentile amplitude so the
    # travelling wavefront is clearly visible at every time.
    vmax = max(np.percentile(np.abs(field), 99.0), 1e-30)
    if have_cartopy:
        ax.set_global()
        ax.gridlines(linewidth=0.5, color="0.7", alpha=0.5)
        ax.pcolormesh(lons, lats, field, transform=ccrs.PlateCarree(),
                      cmap='RdBu_r', vmin=-vmax, vmax=vmax, rasterized=True, shading='auto')
        ax.plot(0, 0, '*', ms=10, color='k', transform=ccrs.PlateCarree(), zorder=10)
    else:
        ax.pcolormesh(lons, lats, field, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      rasterized=True, shading='auto')
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f't = {time[it]:.0f} s', fontsize=14)
fig.suptitle(f'Moon impact: surface vertical displacement ({RES} run, '
             'each panel self-scaled)', fontsize=17)
save_fig(fig, os.path.join(FIG_DIR, 'moon_surface_snapshots_plain'))

# ==== Figure 2: record section (wiggle traces, like ex08 enceladus) ====
print('Moon Fig 2: record section')
fig, ax = plt.subplots(figsize=(7, 8))
# select a subset of receivers spread evenly over distance for clean moveout
n_sel = min(45, len(dist))
sel = np.linspace(0, len(dist) - 1, n_sel).astype(int)
dsel = dist[sel]
span = max(dsel.max() - dsel.min(), 1.0)
gain = span / 30.0           # amplitude of the normalised wiggles, in degrees
for ind, dd in zip(sel, dsel):
    tr = uz[ind]
    nv = max(np.max(np.abs(tr)), 1e-30)
    ax.plot(time, tr / nv * gain + dd, color=COLORS['mono'], lw=0.5)
ax.set_xlabel('Time after impact (s)')
ax.set_ylabel('Epicentral distance (deg)')
ax.set_xlim(time[0], time[-1])
ax.set_ylim(dsel.min() - gain * 1.5, dsel.max() + gain * 1.5)
ax.set_title(f'Moon impact: surface-wave record section ({RES} run)', fontsize=17)
save_fig(fig, os.path.join(FIG_DIR, 'moon_record_section'))

# ==== Figure 3: three-component seismogram ====
print('Moon Fig 3: seismogram')
ncomp = min(3, seis.shape[0])
fig, axes = plt.subplots(ncomp, 1, figsize=(7, 4.5), sharex=True)
axes = np.atleast_1d(axes)
for ich in range(ncomp):
    ch = 'RTZ'[ich]
    axes[ich].plot(seis_time, seis[ich], color=CHANNEL_COLORS[ch], lw=0.7)
    axes[ich].set_ylabel(f'{CHANNEL_LABELS[ch]}\n$U$ (m)', fontsize=15)
    axes[ich].tick_params(direction='in')
axes[0].set_title(f'Moon receiver at Delta = {seis_dist:.0f} deg  ({RES} run)', fontsize=17)
axes[-1].set_xlabel('Time after impact (s)')
axes[0].set_xlim(seis_time[0], seis_time[-1])
save_fig(fig, os.path.join(FIG_DIR, 'moon_seismogram'))

# ==== source-receiver configuration map (Moon globe, no Earth coastlines) ====
# the surface receivers lie on the equator (lat 0) at lon = epicentral distance
# from the on-axis impact at (0, 0).
print('Moon Fig: source-receiver map')
from pub_style import source_receiver_map
_sub = dist[::max(1, len(dist) // 120)]
_recv = [(0.0, float(dd)) for dd in _sub]
source_receiver_map((0.0, 0.0), _recv,
    os.path.join(FIG_DIR, 'moon_source_receiver_map'),
    title='Moon impact — source-receiver geometry', coastlines=False,
    recv_label=f'surface line (n={len(dist)})')

# ==== Combined overview figure: 6 surface snapshots + record section ====
# Arranges the FIRST SIX orthographic surface-displacement snapshots (a 3x2 block
# on the left) beside the wiggle record section (one tall panel on the right) in a
# single balanced figure. Reuses exactly the same data and per-panel logic as the
# stand-alone snapshot and record-section figures above (dist, time, uz).
print('Moon Fig: combined overview')
import matplotlib.gridspec as gridspec

six_t = snap_t[:6]                       # first six snapshot times
n6 = len(six_t)
SR, SC = 3, 2                            # snapshot block: 3 rows x 2 cols

# this figure mixes cartopy GeoAxes with a normal axis; constrained_layout does
# not handle GeoAxes, so disable it up front and place panels manually.
fig = plt.figure(figsize=(11, 8.5))
fig.set_layout_engine('none')
# left half -> 3x2 snapshot grid; right half -> single tall record-section panel.
# width_ratios give the record section a bit more room so its moveout is legible.
outer = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.16, figure=fig)
gs_snap = gridspec.GridSpecFromSubplotSpec(SR, SC, subplot_spec=outer[0],
                                           wspace=0.08, hspace=0.22)

snap_kw = ({'projection': ccrs.Orthographic(0, 0)} if have_cartopy else {})
for k, t in enumerate(six_t):
    ax = fig.add_subplot(gs_snap[k // SC, k % SC], **snap_kw)
    it = int(np.argmin(np.abs(time - t)))
    field = np.interp(DELTA, dist, uz[:, it], left=uz[0, it], right=uz[-1, it])
    vmax = max(np.percentile(np.abs(field), 99.0), 1e-30)
    if have_cartopy:
        ax.set_global()
        ax.gridlines(linewidth=0.5, color="0.7", alpha=0.5)
        ax.pcolormesh(lons, lats, field, transform=ccrs.PlateCarree(),
                      cmap='RdBu_r', vmin=-vmax, vmax=vmax, rasterized=True,
                      shading='auto')
        # epicentral-distance rings from the impact = the spatial scale
        rings = ax.contour(LON, LAT, DELTA, levels=[30, 60, 90], colors='0.2',
                           linewidths=0.5, transform=ccrs.PlateCarree())
        ax.clabel(rings, fmt='%d$^\\circ$', fontsize=9, inline=True, inline_spacing=1)
        ax.plot(0, 0, '*', ms=15, color=COLORS['secondary'], mec='k', mew=0.7,
                transform=ccrs.PlateCarree(), zorder=10)
    else:
        ax.pcolormesh(lons, lats, field, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      rasterized=True, shading='auto')
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f't = {time[it]:.0f} s', fontsize=16)

# right: record section (same selection / gain as the stand-alone Fig 2)
axR = fig.add_subplot(outer[1])
for ind, dd in zip(sel, dsel):
    tr = uz[ind]
    nv = max(np.max(np.abs(tr)), 1e-30)
    axR.plot(time, tr / nv * gain + dd, color=COLORS['mono'], lw=0.5)
axR.set_xlabel('Time after impact (s)', fontsize=19)
axR.set_ylabel('Epicentral distance ($^\\circ$)', fontsize=19)
axR.set_xlim(time[0], time[-1])
axR.set_ylim(dsel.min() - gain * 1.5, dsel.max() + gain * 1.5)
axR.tick_params(labelsize=17)
axR.set_title('Surface-wave record section', fontsize=19)

fig.suptitle('Moon impact — surface wavefield and record section',
             fontsize=22)
fig.subplots_adjust(left=0.04, right=0.97, top=0.91, bottom=0.13)
# wavefield amplitude scale (each snapshot self-scaled to its own peak)
import matplotlib as mpl
cax = fig.add_axes([0.07, 0.065, 0.33, 0.014])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('$U_z$ / max$|U_z|$  (each snapshot normalised)', fontsize=13)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=8)
save_fig(fig, os.path.join(FIG_DIR, 'moon_overview'), tight=False)

print('Moon figures complete.')
