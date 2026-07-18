#!/usr/bin/env python3
"""
Moon impact wavefield on the REAL lunar surface -- 8 time snapshots (still figure).

Same content as moon_wavefield_on_surface.mp4 but as a 2x4 grid of stills: the
bundled lunar imagery (assets/moon_2k.jpg) is the globe texture and the impact
wavefield U(Delta, t) is overlaid, transparent where the ground is quiet so the
wavefront glows over the real surface. Each panel is normalised to its own peak.
"""
import sys, os, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style, save_fig, COLORS
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from PIL import Image

apply_style()
MOON_IMG = os.path.join(HERE, 'assets', 'moon_8k.jpg')
if not os.path.exists(MOON_IMG):
    MOON_IMG = os.path.join(HERE, 'assets', 'moon_2k.jpg')
img = np.asarray(Image.open(MOON_IMG))
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist, time, uz = d['dist'], d['time'], d['uz']

nlat, nlon = 180, 360
lons = np.linspace(-180, 180, nlon); lats = np.linspace(-90, 90, nlat)
LON, LAT = np.meshgrid(lons, lats)
DELTA = np.degrees(np.arccos(np.clip(
    np.cos(np.radians(LAT)) * np.cos(np.radians(LON)), -1, 1)))

snap_t = np.linspace(80, 760, 8)                # 8 times: ring expanding to the limb
NR, NC = 2, 4
fig, axes = plt.subplots(NR, NC, figsize=(15, 8.2),
                         subplot_kw={'projection': ccrs.Orthographic(0, -25)})
axes = axes.ravel()
for ax, ts in zip(axes, snap_t):
    fr = int(np.argmin(np.abs(time - ts)))
    field = np.interp(DELTA, dist, uz[:, fr], left=uz[0, fr], right=uz[-1, fr])
    vmax = max(np.percentile(np.abs(field), 99.0), 1e-30)
    masked = np.ma.array(field, mask=np.abs(field) < 0.18 * vmax)
    ax.set_global()
    ax.imshow(img, origin='upper', extent=[-180, 180, -90, 90],
              transform=ccrs.PlateCarree(), zorder=0)
    ax.pcolormesh(lons, lats, masked, transform=ccrs.PlateCarree(), cmap='RdBu_r',
                  vmin=-vmax, vmax=vmax, alpha=0.78, shading='auto',
                  rasterized=True, zorder=5)
    ax.plot(0, 0, '*', ms=11, color=COLORS['secondary'], mec='k', mew=0.6,
            transform=ccrs.PlateCarree(), zorder=10)
    ax.set_title(f't = {time[fr]:.0f} s', fontsize=16)

cax = fig.add_axes([0.30, 0.07, 0.40, 0.016])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('$U_z$ / max$|U_z|$  (each snapshot normalised)', fontsize=14)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=10)
fig.suptitle('Moon impact — surface wavefield on the lunar surface', fontsize=22)
fig.subplots_adjust(left=0.02, right=0.98, top=0.91, bottom=0.12, wspace=0.05, hspace=0.18)
save_fig(fig, os.path.join(HERE, 'moon_surface_snapshots'), tight=False)
print('Moon real-surface snapshots done.')
