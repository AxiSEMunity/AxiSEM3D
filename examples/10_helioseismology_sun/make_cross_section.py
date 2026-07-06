#!/usr/bin/env python3
"""
Sun helioseismology — combined interior-wavefield overview figure.
  LEFT : the first several meridional cross-section snapshots of the interior
         acoustic wavefield (radial displacement U3), mirrored about the
         rotation axis.
  RIGHT: vertical-displacement record section along the azimuth-0 meridian line
         (epicentral distance vs. time).

This reads the raw AxiSEM3D output that THIS example produces (no pre-reduced
data set):
  output/elements/sun_cross_section/axisem3d_synthetics.nc.rank*
        the meridional element wavefield (U3) with (s, z) coordinates
  output/stations/stationsforplotting/a0.d*.ascii  +  data_time.ascii
        the per-receiver RTZ traces and the shared time axis

The element cross-section group is what makes this the "extreme" run, so this
figure is only produced once you have run the full simulation (see README).
Requires numpy, pandas and netCDF4 in addition to matplotlib.

Run from anywhere:
    python make_cross_section.py
"""
import os, sys, glob
import numpy as np
import pandas as pd
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))

# pub_style.py is the shared plotting-style helper staged in the examples root
# (one level above this example directory).
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import apply_style, save_fig, COLORS
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

apply_style()
R_SUN = 6.96e8                                  # photosphere radius [m]

# INPUT: directory holding this example's AxiSEM3D run output (its ./output/).
INPUT = os.path.join(HERE, 'output')
XS = os.path.join(INPUT, 'elements', 'sun_cross_section')
ST = os.path.join(INPUT, 'stations', 'stationsforplotting')

# ===================== cross-section: first snapshots =====================
fs = sorted(glob.glob(os.path.join(XS, 'axisem3d_synthetics.nc.rank*')))
if not fs:
    sys.exit(
        f"No element cross-section output found in {XS}.\n"
        "Run the simulation first (see README 'Run'); this figure needs the\n"
        "sun_cross_section element group written under output/elements/.")
print(f'{len(fs)} cross-section rank files')
d0 = Dataset(fs[0]); d0.set_auto_mask(False)
t = np.asarray(d0.variables['data_time'][:]); d0.close()
NS = 6
snap_idx = np.linspace(8, len(t) - 1, 8).astype(int)[:NS]    # the FIRST six
print('snapshot times (h):', [f'{t[i]/3600:.1f}' for i in snap_idx])

S, Z, W = [], [], []
GC = 4                                          # centre GLL of the 3x3 element
for k, f in enumerate(fs):
    d = Dataset(f); d.set_auto_mask(False)
    co = np.asarray(d.variables['list_element_coords'][:])           # (elem,9,2)
    w = np.asarray(d.variables['data_wave__NaG=1'][:])               # (elem,1,9,1,time)
    S.append(co[:, GC, 0]); Z.append(co[:, GC, 1])
    W.append(w[:, 0, GC, 0, :][:, snap_idx])
    d.close()
S = np.concatenate(S); Z = np.concatenate(Z); W = np.concatenate(W, axis=0)
Sfull = np.concatenate([S, -S]) / 1e8           # mirror about the axis; units 1e8 m
Zfull = np.concatenate([Z, Z]) / 1e8
print(f'{S.size:,} element centres')

# ===================== record section: a0 meridian line =====================
time_st = np.loadtxt(os.path.join(ST, 'data_time.ascii'))
DS = max(1, len(time_st) // 3000)               # downsample for a light wiggle plot
ts = time_st[::DS]
recs = []
for f in sorted(glob.glob(os.path.join(ST, 'a0.d*.ascii'))):
    dist = int(os.path.basename(f).split('.d')[1].split('.ascii')[0])
    z = pd.read_csv(f, sep=r'\s+', header=None, usecols=[2]).values[::DS, 0]
    recs.append((dist, z))
recs.sort()
dists = np.array([r[0] for r in recs])
gain = (dists.max() - dists.min()) / 30.0       # wiggle amplitude in degrees
print(f'{len(recs)} a0 receivers, dist {dists.min()}-{dists.max()} deg')

# ===================== figure: snapshots (left) + record section (right) ====
fig = plt.figure(figsize=(11, 8.5))
fig.set_layout_engine('none')
outer = gridspec.GridSpec(1, 2, width_ratios=[1.15, 1.0], wspace=0.14, figure=fig)
gs_snap = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=outer[0],
                                           wspace=0.05, hspace=0.18)
theta = np.linspace(0, 2 * np.pi, 400)
for j in range(NS):
    ax = fig.add_subplot(gs_snap[j // 2, j % 2])
    val = W[:, j]
    vmax = max(np.percentile(np.abs(val), 98.0), 1e-30)
    ax.scatter(Sfull, Zfull, c=np.concatenate([val, val]), s=0.6, cmap='RdBu_r',
               vmin=-vmax, vmax=vmax, rasterized=True, linewidths=0, marker='s')
    ax.plot(R_SUN / 1e8 * np.cos(theta), R_SUN / 1e8 * np.sin(theta),
            color='#555555', lw=0.7)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f't = {t[snap_idx[j]] / 3600:.1f} h', fontsize=16)
    if j == 0:                                  # 200 Mm scale bar on the first panel
        ax.plot([-1.0, 1.0], [-5.6, -5.6], color='k', lw=2.0, solid_capstyle='butt')
        ax.text(0, -4.9, '200 Mm', ha='center', va='bottom', fontsize=8)

# right: vertical-displacement record section
axR = fig.add_subplot(outer[1])
for dist, z in recs:
    nv = max(np.max(np.abs(z)), 1e-30)
    axR.plot(ts, z / nv * gain + dist, color=COLORS['mono'], lw=0.6)
axR.set_xlabel('Time after source origin (s)', fontsize=17)
axR.set_ylabel('Epicentral distance (deg)', fontsize=17)
axR.set_xlim(ts[0], ts[-1])
axR.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
axR.tick_params(labelsize=15)
axR.set_title('Vertical-displacement record section', fontsize=17)

# amplitude scale for the snapshots (diverging, each panel normalised to its peak)
import matplotlib as mpl
cax = fig.add_axes([0.07, 0.045, 0.34, 0.014])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label(r'$U_3\,/\,\max|U_3|$  (each snapshot normalised)', fontsize=13)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=8)

fig.suptitle('Sun helioseismology - interior wavefield cross-section and record section',
             fontsize=20)
fig.subplots_adjust(left=0.02, right=0.97, top=0.92, bottom=0.13)
save_fig(fig, os.path.join(HERE, 'figures', 'sun_overview'), tight=False)
print('Sun overview done.')
