#!/usr/bin/env python3
"""
Moon impact wavefield draped on the REAL lunar surface (orthographic globe).

Uses the bundled lunar imagery (assets/moon_2k.jpg -- the actual albedo map: maria,
highlands, craters) as the globe texture, and overlays the axisymmetric impact
wavefield U(Delta, t) from ./data/moon_surface.npz, made transparent where the
ground is quiet so the wavefront glows over the real surface. Writes one frame; set
MOVIE=True for a GIF (Pillow, ffmpeg-free).
"""
import sys, os, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style, COLORS
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
from PIL import Image

apply_style()
MOVIE = '--movie' in sys.argv
MOON_IMG = os.path.join(HERE, 'assets', 'moon_8k.jpg')
if not os.path.exists(MOON_IMG):
    MOON_IMG = os.path.join(HERE, 'assets', 'moon_2k.jpg')
img = np.asarray(Image.open(MOON_IMG))                  # equirectangular albedo
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist, time, uz = d['dist'], d['time'], d['uz']

nlat, nlon = 180, 360
lats = np.linspace(-90, 90, nlat); lons = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lons, lats)
DELTA = np.degrees(np.arccos(np.clip(
    np.cos(np.radians(LAT)) * np.cos(np.radians(LON)), -1, 1)))


def field_at(fr):
    return np.interp(DELTA, dist, uz[:, fr], left=uz[0, fr], right=uz[-1, fr])


fig = plt.figure(figsize=(5.2, 5.6))
ax = fig.add_subplot(projection=ccrs.Orthographic(0, -25))   # tilt to see the ring
ax.set_global()
# the real Moon surface as the base layer
ax.imshow(img, origin='upper', extent=[-180, 180, -90, 90],
          transform=ccrs.PlateCarree(), zorder=0)
ax.plot(0, 0, '*', ms=14, color=COLORS['secondary'], mec='k', mew=0.7,
        transform=ccrs.PlateCarree(), zorder=10)
ttl = ax.set_title('', fontsize=14)
state = {'m': None}


def draw(fr):
    field = field_at(fr)
    vmax = max(np.percentile(np.abs(field), 99.0), 1e-30)
    # transparent where the ground is quiet -> the wavefront glows over the surface
    masked = np.ma.array(field, mask=np.abs(field) < 0.18 * vmax)
    if state['m'] is not None:
        state['m'].remove()
    state['m'] = ax.pcolormesh(lons, lats, masked, transform=ccrs.PlateCarree(),
                               cmap='RdBu_r', vmin=-vmax, vmax=vmax, alpha=0.75,
                               shading='auto', rasterized=True, zorder=5)
    ttl.set_text(f'Moon impact wavefield on the lunar surface   t = {time[fr]:.0f} s')


if MOVIE:
    import imageio_ffmpeg
    plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
    frames = list(range(0, len(time), 6))
    ani = animation.FuncAnimation(fig, draw, frames=frames, blit=False)
    mdir = os.path.join(HERE, 'movies'); os.makedirs(mdir, exist_ok=True)
    out = os.path.join(mdir, 'moon_wavefield_on_surface.mp4')
    ani.save(out, writer=animation.FFMpegWriter(fps=14, bitrate=6000), dpi=120)
    print('Saved', out)
else:
    draw(int(np.argmin(np.abs(time - 330))))            # a nicely-expanded ring
    out = os.path.join(HERE, 'moon_wavefield_on_surface.png')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print('Saved', out)
