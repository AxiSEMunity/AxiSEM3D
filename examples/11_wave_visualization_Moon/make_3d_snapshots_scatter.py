#!/usr/bin/env python3
"""
Moon impact WITH scattering -- wavefield on the 3-D textured lunar globe, 8 snapshots.

Unlike the smooth (axisymmetric) run, the scattered field is genuinely 2-D on the
surface, so it is read from a lat/lon grid of receivers (surface_grid station group,
U3) rather than a single meridian line. Combines the per-rank NETCDF station files
into U(lat, lon, t), interpolates each snapshot onto the PyVista globe, and drapes it
over the real lunar texture -- the direct analogue of make_3d_snapshots.py but for the
scattering run, so the two figures can be compared side by side.

Requires the megaregolith-scattering run's own output, which is NOT bundled here
(too large): run the solver on sim/input_scatter/ (see FIGURES.md) so that
sim/output_scatter/stations/surface_grid/ exists, then run this script.
"""
import os, sys, glob
import numpy as np
import pyvista as pv
from matplotlib import cm
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

pv.OFF_SCREEN = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style
apply_style()
mpl.rcParams.update({'pdf.use14corefonts': True, 'font.family': 'Helvetica',
                     'axes.unicode_minus': False})

MOON_IMG = os.path.join(HERE, 'assets', 'moon_8k.jpg')
if not os.path.exists(MOON_IMG):
    MOON_IMG = os.path.join(HERE, 'assets', 'moon_2k.jpg')
GRID_TXT = os.path.join(HERE, 'sim', 'input_scatter', 'moon_grid.txt')
STA_DIR = os.path.join(HERE, 'sim', 'output_scatter', 'stations', 'surface_grid')

if not os.path.exists(GRID_TXT):
    raise SystemExit(
        f'{GRID_TXT} not found -- this script needs the megaregolith-scattering '
        'run\'s input/output (see FIGURES.md for how to obtain/reproduce it).')

# ---- receiver grid: station name -> (lat, lon) ----
g = np.loadtxt(GRID_TXT, dtype=str)
gname = g[:, 0]; glat = g[:, 2].astype(float); glon = g[:, 3].astype(float)
lats_g = np.unique(glat); lons_g = np.unique(glon)
row = {n: (int(np.searchsorted(lats_g, glat[i])), int(np.searchsorted(lons_g, glon[i])))
       for i, n in enumerate(gname)}


def gkey(raw):
    s = raw.strip()
    return s.split('.')[-1] if '.' in s else s          # 'GR.G00042' -> 'G00042'


# ---- read + combine per-rank NETCDF station files -> Ug(lat, lon, time) ----
fs = sorted(glob.glob(os.path.join(STA_DIR, 'axisem3d_synthetics.nc.rank*')))
if not fs:
    raise SystemExit('no surface_grid output yet -- run must finish first')
d0 = Dataset(fs[0]); d0.set_auto_mask(False)
t = np.asarray(d0.variables['data_time'][:]); d0.close()
Ug = np.full((len(lats_g), len(lons_g), len(t)), np.nan, dtype='f4')
for f in fs:
    d = Dataset(f); d.set_auto_mask(False)
    w = np.asarray(d.variables['data_wave'][:])          # (nsta, nchan, ntime)
    names = [gkey(b''.join(r).decode('latin1')) for r in np.asarray(d.variables['list_station'][:])]
    for i, n in enumerate(names):
        if n in row:
            r, c = row[n]; Ug[r, c] = w[i, 0]            # channel 0 = U3 (vertical)
    d.close()
print(f'{np.isfinite(Ug[:, :, 0]).sum()} grid stations, {len(t)} frames, '
      f't = {t[0]:.0f}..{t[-1]:.0f} s')
# the grid cell directly under the source carries the static near-field
# offset (~1e5-1e6x the surrounding wavefield) -- bilinear interpolation
# smears that single spike into a saturated square. Replace it with its
# 4 neighbours so only the propagating wavefield is drawn.
src_r = int(np.argmin(np.abs(lats_g)))
src_c = int(np.argmin(np.abs(lons_g)))
nbrs = [(src_r + dr, src_c + dc) for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]
        if 0 <= src_r + dr < len(lats_g)]
Ug[src_r, src_c] = np.nanmean([Ug[r, c % len(lons_g)] for r, c in nbrs], axis=0)
# fill any gaps + wrap longitude for interpolation
Ug = np.nan_to_num(Ug)
lons_ext = np.r_[lons_g, lons_g[0] + 360]
Ug_ext = np.concatenate([Ug, Ug[:, :1]], axis=1)

# ---- PyVista textured globe (as in make_3d_snapshots) ----
R, nlat, nlon = 1.0, 481, 961
la = np.linspace(-90, 90, nlat); lo = np.linspace(-180, 180, nlon)
LO, LA = np.meshgrid(lo, la)
lr, gr = np.radians(LA), np.radians(LO)
X = np.cos(lr) * np.cos(gr); Y = np.cos(lr) * np.sin(gr); Z = np.sin(lr)
base = pv.StructuredGrid(X, Y, Z)
bp = base.points
plon = np.degrees(np.arctan2(bp[:, 1], bp[:, 0]))
plat = np.degrees(np.arcsin(np.clip(bp[:, 2], -1, 1)))
base.active_texture_coordinates = np.c_[(plon + 180) / 360.0, (plat + 90) / 180.0]
tex = pv.Texture(MOON_IMG)
ov = pv.StructuredGrid(1.004 * X, 1.004 * Y, 1.004 * Z)
qlon = np.where(plon < lons_g[0], plon + 360, plon)      # for the wrapped grid
qpts = np.c_[np.clip(plat, lats_g[0], lats_g[-1]), np.clip(qlon, lons_ext[0], lons_ext[-1])]


def render(fr):
    rgi = RegularGridInterpolator((lats_g, lons_ext), Ug_ext[:, :, fr],
                                  bounds_error=False, fill_value=0.0)
    U = rgi(qpts)
    vmax = max(np.percentile(np.abs(U), 99.0), 1e-30)
    rgba = cm.RdBu_r(Normalize(-vmax, vmax)(U))
    rgba[:, 3] = np.clip(np.abs(U) / vmax, 0, 1) ** 0.7
    ov['rgba'] = (rgba * 255).astype(np.uint8)
    pl = pv.Plotter(off_screen=True, window_size=(1000, 1000)); pl.set_background('white')
    pl.add_mesh(base, texture=tex, smooth_shading=True, ambient=0.35, diffuse=0.8, specular=0.05)
    pl.add_mesh(ov, scalars='rgba', rgba=True, show_scalar_bar=False,
                smooth_shading=True, reset_camera=False)
    pl.camera_position = [(4.0, -1.7, 1.3), (0, 0, 0), (0, 0, 1)]
    img = np.asarray(pl.screenshot(return_img=True)); pl.close()
    return img


snap_t = np.linspace(80, 900, 8)
shots = [int(np.argmin(np.abs(t - ts))) for ts in snap_t]
imgs = [render(fr) for fr in shots]

mpl.rcParams['figure.constrained_layout.use'] = False
fig, axes = plt.subplots(2, 4, figsize=(16, 8.6))
for ax, im, fr in zip(axes.ravel(), imgs, shots):
    ax.imshow(im); ax.axis('off'); ax.set_title(f't = {t[fr]:.0f} s', fontsize=17)
cax = fig.add_axes([0.32, 0.06, 0.36, 0.016])
cb = fig.colorbar(cm.ScalarMappable(norm=Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('Uz / max|Uz|  (each snapshot normalised)', fontsize=14)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=10)
fig.suptitle('Moon impact with scattering — wavefield on the 3-D lunar surface', fontsize=22, y=0.995)
fig.subplots_adjust(left=0.01, right=0.99, top=0.85, bottom=0.11, wspace=0.02, hspace=0.25)
fig.savefig(os.path.join(HERE, 'moon_3d_snapshots_scatter.pdf'), dpi=200)
fig.savefig(os.path.join(HERE, 'moon_3d_snapshots_scatter.png'), dpi=160)
print('Saved moon_3d_snapshots_scatter.pdf / .png')
