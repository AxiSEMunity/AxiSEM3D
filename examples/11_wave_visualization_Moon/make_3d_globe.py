#!/usr/bin/env python3
"""
Moon impact wavefield on a TRUE 3-D textured lunar globe (PyVista, headless).

Base sphere is textured with the real lunar imagery (assets/moon_8k.jpg). A thin
overlay shell (radius 1.004) carries the axisymmetric impact field U(Delta, t),
coloured RdBu_r and made transparent where the ground is quiet (nan_opacity=0) so
the wavefront glows over the real surface in 3-D. Renders a still by default;
pass --spin for a rotating MP4, --time for a time-evolving MP4.
"""
import os, sys, glob
import numpy as np
import pyvista as pv
from matplotlib import cm
from matplotlib.colors import Normalize

pv.OFF_SCREEN = True
HERE = os.path.dirname(os.path.abspath(__file__))
MOON_IMG = os.path.join(HERE, 'assets', 'moon_8k.jpg')
if not os.path.exists(MOON_IMG):
    MOON_IMG = os.path.join(HERE, 'assets', 'moon_2k.jpg')
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist, time, uz = d['dist'], d['time'], d['uz']

# ---- base textured sphere (real Moon) ----
R, nlat, nlon = 1.0, 481, 961
lat = np.linspace(-90, 90, nlat); lon = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lon, lat)
lr, gr = np.radians(LAT), np.radians(LON)
X = R * np.cos(lr) * np.cos(gr); Y = R * np.cos(lr) * np.sin(gr); Z = R * np.sin(lr)
base = pv.StructuredGrid(X, Y, Z)
p = base.points
plon = np.degrees(np.arctan2(p[:, 1], p[:, 0]))
plat = np.degrees(np.arcsin(np.clip(p[:, 2] / R, -1, 1)))
base.active_texture_coordinates = np.c_[(plon + 180) / 360.0, (plat + 90) / 180.0]
tex = pv.Texture(MOON_IMG)

# ---- overlay shell carrying the wavefield ----
Ro = 1.004
ov = pv.StructuredGrid(Ro * X, Ro * Y, Ro * Z)
ovdir = ov.points / Ro
IMPACT = np.array([1.0, 0.0, 0.0])               # impact faces +x (toward camera)
delta = np.degrees(np.arccos(np.clip(ovdir @ IMPACT, -1, 1)))


def wavefield(fr):
    U = np.interp(delta, dist, uz[:, fr], left=uz[0, fr], right=uz[-1, fr])
    vmax = max(np.percentile(np.abs(U), 99.0), 1e-30)
    rgba = cm.RdBu_r(Normalize(-vmax, vmax)(U))          # per-point colour
    rgba[:, 3] = np.clip(np.abs(U) / vmax, 0, 1) ** 0.7  # alpha = wavefront fade
    return (rgba * 255).astype(np.uint8), vmax


def make_plotter():
    pl = pv.Plotter(off_screen=True, window_size=(1500, 1500))
    pl.set_background('white')
    pl.add_mesh(base, texture=tex, smooth_shading=True, ambient=0.35,
                diffuse=0.8, specular=0.05)
    pl.camera_position = [(4.0, -1.7, 1.3), (0, 0, 0), (0, 0, 1)]
    return pl


MODE = 'spin' if '--spin' in sys.argv else ('time' if '--time' in sys.argv else 'still')

if MODE == 'still':
    fr = int(np.argmin(np.abs(time - 330)))
    rgba, vmax = wavefield(fr)
    ov['rgba'] = rgba
    pl = make_plotter()
    pl.add_mesh(ov, scalars='rgba', rgba=True, show_scalar_bar=False,
                smooth_shading=True, reset_camera=False)
    pl.add_text(f'Moon impact on the lunar surface   t = {time[fr]:.0f} s',
                position='upper_edge', font_size=14, color='black')
    out = os.path.join(HERE, 'moon_3d_globe.png')
    pl.screenshot(out, transparent_background=False)
    print('Saved', out)
else:
    import imageio_ffmpeg
    mdir = os.path.join(HERE, 'movies'); os.makedirs(mdir, exist_ok=True)
    out = os.path.join(mdir, 'moon_3d_%s.mp4' % MODE)
    wr = imageio_ffmpeg.write_frames(out, (1500, 1500), fps=20, quality=7)
    wr.send(None)
    if MODE == 'spin':
        fr = int(np.argmin(np.abs(time - 330))); rgba, vmax = wavefield(fr)
        ov['rgba'] = rgba
        pl = make_plotter()
        pl.add_mesh(ov, scalars='rgba', rgba=True, show_scalar_bar=False,
                    reset_camera=False)
        for az in range(0, 360, 3):
            pl.camera.azimuth = az
            wr.send(np.asarray(pl.screenshot(return_img=True))[:, :, :3].copy())
    else:                                            # time evolution
        pl = make_plotter(); actor = None
        for fr in range(0, len(time), 6):
            rgba, vmax = wavefield(fr); ov['rgba'] = rgba
            if actor is not None:
                pl.remove_actor(actor)
            actor = pl.add_mesh(ov, scalars='rgba', rgba=True,
                                show_scalar_bar=False, reset_camera=False)
            wr.send(np.asarray(pl.screenshot(return_img=True))[:, :, :3].copy())
    wr.close()
    print('Saved', out)
