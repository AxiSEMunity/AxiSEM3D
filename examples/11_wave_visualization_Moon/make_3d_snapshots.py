#!/usr/bin/env python3
"""
Moon impact wavefield on the TRUE 3-D textured lunar globe -- 8 time snapshots.

Renders 8 PyVista globe views (real lunar texture + wavefield RGBA overlay, the same
look as moon_3d_globe.png) at successive times and composites them into one 2x4
publication figure (Helvetica, shared colorbar). The 3-D analogue of
moon_surface_snapshots_real.pdf.
"""
import os, sys, glob
import numpy as np
import pyvista as pv
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib as mpl

pv.OFF_SCREEN = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style
apply_style()
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['pdf.use14corefonts'] = True    # reference core Helvetica (Illustrator-editable)
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['axes.unicode_minus'] = False

MOON_IMG = glob.glob(os.path.join(HERE, '..', '..',
                     '09_wave_visualization_Moon', 'assets', 'moon_8k.jpg'))[0]
d = np.load(os.path.join(HERE, 'data', 'moon_surface.npz'))
dist, time, uz = d['dist'], d['time'], d['uz']

# ---- base textured sphere + wavefield overlay shell (as in make_3d_globe) ----
R, nlat, nlon = 1.0, 481, 961
lat = np.linspace(-90, 90, nlat); lon = np.linspace(-180, 180, nlon)
LON, LAT = np.meshgrid(lon, lat)
lr, gr = np.radians(LAT), np.radians(LON)
X = R * np.cos(lr) * np.cos(gr); Y = R * np.cos(lr) * np.sin(gr); Z = R * np.sin(lr)
base = pv.StructuredGrid(X, Y, Z)
bp = base.points
plon = np.degrees(np.arctan2(bp[:, 1], bp[:, 0]))
plat = np.degrees(np.arcsin(np.clip(bp[:, 2] / R, -1, 1)))
base.active_texture_coordinates = np.c_[(plon + 180) / 360.0, (plat + 90) / 180.0]
tex = pv.Texture(MOON_IMG)
Ro = 1.004
ov = pv.StructuredGrid(Ro * X, Ro * Y, Ro * Z)
delta = np.degrees(np.arccos(np.clip((ov.points / Ro) @ np.array([1., 0., 0.]), -1, 1)))


def render(fr):
    U = np.interp(delta, dist, uz[:, fr], left=uz[0, fr], right=uz[-1, fr])
    vmax = max(np.percentile(np.abs(U), 99.0), 1e-30)
    rgba = cm.RdBu_r(Normalize(-vmax, vmax)(U))
    rgba[:, 3] = np.clip(np.abs(U) / vmax, 0, 1) ** 0.7
    ov['rgba'] = (rgba * 255).astype(np.uint8)
    pl = pv.Plotter(off_screen=True, window_size=(1000, 1000))
    pl.set_background('white')
    pl.add_mesh(base, texture=tex, smooth_shading=True, ambient=0.35,
                diffuse=0.8, specular=0.05)
    pl.add_mesh(ov, scalars='rgba', rgba=True, show_scalar_bar=False,
                smooth_shading=True, reset_camera=False)
    pl.camera_position = [(4.0, -1.7, 1.3), (0, 0, 0), (0, 0, 1)]
    img = pl.screenshot(return_img=True)
    pl.close()
    return np.asarray(img)


snap_t = np.linspace(80, 760, 8)               # ring expanding to the limb
print('rendering 8 globes...')
shots = [(int(np.argmin(np.abs(time - ts))),) for ts in snap_t]
imgs = [render(fr) for (fr,) in shots]

# ---- composite into a 2x4 figure ----
mpl.rcParams['figure.constrained_layout.use'] = False
fig, axes = plt.subplots(2, 4, figsize=(16, 8.6))
for ax, img, (fr,) in zip(axes.ravel(), imgs, shots):
    ax.imshow(img); ax.axis('off')
    ax.set_title(f't = {time[fr]:.0f} s', fontsize=17)

cax = fig.add_axes([0.32, 0.06, 0.36, 0.016])
cb = fig.colorbar(cm.ScalarMappable(norm=Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('Uz / max|Uz|  (each snapshot normalised)', fontsize=14)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=10)
fig.suptitle('Moon impact — wavefield on the 3-D lunar surface', fontsize=22, y=0.995)
fig.subplots_adjust(left=0.01, right=0.99, top=0.85, bottom=0.11, wspace=0.02, hspace=0.25)
fig.savefig(os.path.join(HERE, 'moon_3d_snapshots.pdf'), dpi=200)
fig.savefig(os.path.join(HERE, 'moon_3d_snapshots.png'), dpi=160)
print('Saved moon_3d_snapshots.pdf / .png')
