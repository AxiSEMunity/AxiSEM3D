#!/usr/bin/env python3
"""
Moon meteoroid-impact SURFACE wavefield -- animation (orthographic lunar globe).

Reuses exactly the per-frame reconstruction of the still figures: the axisymmetric
impact field U(Delta, t) along one meridian (./data/moon_surface.npz) is mapped onto
every (lat, lon) by its epicentral distance Delta from the impact. Each frame is
normalised to its own peak so the wavefront stays visible as it spreads and decays.
Writes moon_wavefield.gif (ffmpeg-free, via Pillow).
"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style, COLORS
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs

apply_style()
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist, time, uz = d['dist'], d['time'], d['uz']
RES = str(d['res'])

nlat, nlon = 140, 280
lats = np.linspace(-90, 90, nlat); lons = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lons, lats)
DELTA = np.degrees(np.arccos(np.clip(
    np.cos(np.radians(LAT)) * np.cos(np.radians(LON)), -1, 1)))

frames = list(range(0, len(time), 5))           # 900 -> ~180 frames
print(f'{len(frames)} frames, t = {time[0]:.0f}..{time[-1]:.0f} s ({RES} run)')


def field_at(fr):
    return np.interp(DELTA, dist, uz[:, fr], left=uz[0, fr], right=uz[-1, fr])


fig = plt.figure(figsize=(5.0, 5.4))
ax = fig.add_subplot(projection=ccrs.Orthographic(0, 0))
ax.set_global()
ax.gridlines(linewidth=0.4, color='0.7', alpha=0.5)
ax.plot(0, 0, '*', ms=13, color=COLORS['secondary'], mec='k', mew=0.6,
        transform=ccrs.PlateCarree(), zorder=10)
ttl = ax.set_title('', fontsize=15)
state = {'mesh': None}


def update(fr):
    field = field_at(fr)
    vmax = max(np.percentile(np.abs(field), 99.0), 1e-30)
    if state['mesh'] is not None:
        state['mesh'].remove()
    state['mesh'] = ax.pcolormesh(lons, lats, field, transform=ccrs.PlateCarree(),
                                  cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                                  shading='auto', rasterized=True, zorder=1)
    ttl.set_text(f'Moon impact — surface vertical displacement   t = {time[fr]:.0f} s')
    return state['mesh'], ttl


import imageio_ffmpeg
plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
ani = animation.FuncAnimation(fig, update, frames=frames, blit=False)
mdir = os.path.join(HERE, 'movies'); os.makedirs(mdir, exist_ok=True)
out = os.path.join(mdir, 'moon_wavefield.mp4')
ani.save(out, writer=animation.FFMpegWriter(fps=15, bitrate=6000), dpi=110)
print('Saved', out)
