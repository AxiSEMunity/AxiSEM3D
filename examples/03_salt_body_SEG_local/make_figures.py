#!/usr/bin/env python3
"""
Example 03: Local SEG/EAGE C3 salt body (5 Hz) -- publication figures.

  Fig 1: ex03_ocean_floor_waveforms     -- 3-component (R/T/Z) waveforms at
                                            three ocean-floor locations.
  Fig 2: ex03_cross_section_snapshots    -- 2x4 time-lapse of the vertical
                                            displacement U_Z on the phi = 0
                                            in-plane cross-section.
  Fig 3: ex03_salt_model_velocities      -- the INPUT SEG/EAGE C3 velocity
                                            model (optional; needs an extract
                                            of the large input model, see the
                                            note next to SALT_MODEL_NPZ below).

All inputs are located RELATIVE to this script, so it runs unchanged after you
clone the repository and run the example. The shared plotting style lives in
examples/pub_style.py (one directory up); we add that directory to sys.path
relatively below.
"""
import sys
import os

import numpy as np
import xarray as xr

# --- locate the example dir and the shared style helper, both relatively ---
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # examples/  -> pub_style.py
from pub_style import apply_style, save_fig, CHANNEL_COLORS, CHANNEL_LABELS
import matplotlib.pyplot as plt                     # noqa: E402  (Agg set in pub_style)
import matplotlib.tri as mtri                       # noqa: E402

apply_style()

# ---------------------------------------------------------------------------
# INPUT: the AxiSEM3D output produced by running the 3D simulation in this
# example (see the README "Reproduce" section). Override by setting
# the EX03_OUTPUT environment variable, otherwise default to ./simu3D/output.
# ---------------------------------------------------------------------------
OUTPUT_DIR = os.environ.get('EX03_OUTPUT', os.path.join(HERE, 'simu3D', 'output'))
ELEMENTS_DIR = os.path.join(OUTPUT_DIR, 'elements')
FIG_DIR = os.path.join(HERE, 'figures')

# Optional input for Fig 3: a compact (~0.4 MB) extract of two slices (one
# vertical, one horizontal) of the full ~1 GB SEG_C3_SOLID.nc input model.
# It is NOT produced by simply running the example -- see the README figure
# section for how to derive it. Fig 3 is skipped if this file is absent.
SALT_MODEL_NPZ = os.path.join(HERE, 'data', 'salt_model.npz')


def read_element_output(data_dir):
    """Read an AxiSEM3D element-wise NetCDF group (one file per MPI rank)."""
    nc_fnames = sorted([f for f in os.listdir(data_dir)
                        if 'axisem3d_synthetics.nc' in f])
    print(f'  Reading {len(nc_fnames)} file(s) from {data_dir}')
    nc_files = [xr.open_dataset(os.path.join(data_dir, fn)) for fn in nc_fnames]
    na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)
    data_time = nc_files[0].data_vars['data_time'].values

    xda_na = xr.concat([f.data_vars['list_element_na'] for f in nc_files],
                       dim='dim_element')
    xda_coords = xr.concat([f.data_vars['list_element_coords'] for f in nc_files],
                           dim='dim_element')
    dict_data_wave = {}
    for nag in na_grid:
        dict_data_wave[nag] = xr.concat(
            [f.data_vars[f'data_wave__NaG={nag}'] for f in nc_files],
            dim=f'dim_element__NaG={nag}').values

    return na_grid, data_time, xda_na.values.astype(int), xda_coords.values, dict_data_wave


def get_waveform_xy(xy, na_grid, list_element_na, list_element_coords, dict_data_wave):
    """Reconstruct the time series at horizontal location (x, y) from the
    stored azimuthal Fourier coefficients (inplane + Fourier interpolation)."""
    s = np.linalg.norm(xy)
    phi = np.arctan2(xy[1], xy[0])
    if s > np.max(list_element_coords[:, :, 0]):
        return None
    s_center = list_element_coords[:, 2, 0]
    idx = np.argmin(np.abs(s - s_center))
    s_elem = list_element_coords[idx, :, 0]
    iA = np.argmin(np.abs(s - s_elem))
    s_copy = s_elem.copy()
    s_copy[iA] = 1e100
    iB = np.argmin(np.abs(s - s_copy))
    fA = (s_elem[iB] - s) / (s_elem[iB] - s_elem[iA])
    fB = 1 - fA
    ena = list_element_na[idx]
    dA = dict_data_wave[ena[2]][ena[4], :, iA]
    dB = dict_data_wave[ena[2]][ena[4], :, iB]
    data_wave = dA * fA + dB * fB
    max_order = ena[1] // 2
    result = data_wave[0].copy()
    for order in range(1, max_order + 1):
        coeff = np.zeros(result.shape, dtype=np.complex128)
        coeff.real = data_wave[order * 2 - 1]
        if order * 2 < len(data_wave):
            coeff.imag = data_wave[order * 2]
        result += (2.0 * np.exp(1j * order * phi) * coeff).real
    return result


# ==== Figure 1: Ocean floor waveforms ====
print('Ex03 Fig 1: Ocean floor waveforms')
data_dir = os.path.join(ELEMENTS_DIR, 'Fourier_coefficients_ocean_floor')
na_grid, data_time, el_na, el_coords, data_wave = read_element_output(data_dir)

fig, axes = plt.subplots(3, 1, figsize=(7, 4.5), sharex=True, sharey=True)
for xy, ls, lbl in [([2000, 3000], '-', '(2, 3) km'),
                    ([5000, 0], '--', '(5, 0) km'),
                    ([0, 5000], ':', '(0, 5) km')]:
    w = get_waveform_xy(xy, na_grid, el_na, el_coords, data_wave)
    if w is None:
        continue
    for ich, ch in enumerate('RTZ'):
        axes[ich].plot(data_time, w[ich], ls=ls, lw=0.8,
                       color=CHANNEL_COLORS[ch], label=lbl if ich == 0 else None)

for ich, ch in enumerate('RTZ'):
    axes[ich].set_ylabel(f'{CHANNEL_LABELS[ch]} (m)', fontsize=9)
    axes[ich].tick_params(direction='in')
    if ich < 2:
        axes[ich].tick_params(labelbottom=False)
axes[0].legend(loc='upper right', fontsize=8, title='Location', title_fontsize=8)
axes[0].set_title('SEG Salt Body -- Ocean floor waveforms', fontsize=11)
axes[-1].set_xlabel('Time (s)')
axes[0].set_xlim(0, data_time[-1])  # source-time function is acausal; start at origin
save_fig(fig, os.path.join(FIG_DIR, 'ex03_ocean_floor_waveforms'))

# ==== Figure 2: Inplane cross-section snapshots ====
print('Ex03 Fig 2: Inplane cross-sections')
data_dir2 = os.path.join(ELEMENTS_DIR, 'orthogonal_azimuthal_slices')
na_grid2, data_time2, el_na2, el_coords2, data_wave2 = read_element_output(data_dir2)

nelem = el_coords2.shape[0]
ngll = el_coords2.shape[1]
ntime = len(data_time2)
phi_slices = np.radians(np.arange(0, 360, 90))
nslice = len(phi_slices)

coords_sz = el_coords2.reshape((nelem * ngll, 2))

wave_dim = 2   # Z component of the RTZ frame
islice = 0
phi = phi_slices[islice]

wave = np.ndarray((nelem * ngll, ntime))
for ie in range(nelem):
    wave[ie * ngll:ie * ngll + ngll, :] = data_wave2[nslice][ie, islice, :, wave_dim, :]

s_coords = coords_sz[:, 0] * np.cos(phi)
z_coords = coords_sz[:, 1]
# convert spherical radius -> depth below the top of the imaged model (m)
depth_coords = z_coords.max() - z_coords

# triangulate the GLL points once (geometry is static across snapshots)
triang = mtri.Triangulation(s_coords, depth_coords)

NR, NC = 2, 4   # 8-panel time-lapse
snapshots = [int(ntime * f) for f in np.linspace(0.08, 0.85, NR * NC)]
snapshots = [t for t in snapshots if t < ntime]

vmax = np.percentile(np.abs(wave), 99.5)
levels = np.linspace(-vmax, vmax, 41)

fig, axes = plt.subplots(NR, NC, figsize=(3.0 * NC, 3.3 * NR),
                         sharex=True, sharey=True, constrained_layout=True)
axes = np.atleast_1d(axes).ravel()
for i, tstep in enumerate(snapshots):
    ax = axes[i]
    cf = ax.tricontourf(triang, wave[:, tstep], levels=levels,
                        cmap='RdBu_r', extend='both')
    cf.set_rasterized(True)  # keep the PDF small
    ax.set_title(f't = {data_time2[tstep]:.2f} s', fontsize=9)
    ax.set_aspect('equal')
    ax.set_xlim(s_coords.min(), s_coords.max())
    ax.set_ylim(depth_coords.max(), depth_coords.min())  # depth increases downward
    if i % NC == 0:
        ax.set_ylabel('Depth below top (m)', fontsize=9)
    if i >= len(snapshots) - NC:
        ax.set_xlabel('Distance (m)', fontsize=9)
for k in range(len(snapshots), len(axes)):
    axes[k].set_visible(False)

cbar = fig.colorbar(cf, ax=axes.tolist(), shrink=0.6, pad=0.02)
cbar.set_label('$U_Z$ (m)', fontsize=9)
fig.suptitle('SEG Salt Body -- Vertical displacement cross-section time-lapse',
             fontsize=12)
save_fig(fig, os.path.join(FIG_DIR, 'ex03_cross_section_snapshots'), tight=False)

# ==== Figure 3: the SEG/EAGE C3 salt body INPUT velocity model (optional) ====
if os.path.isfile(SALT_MODEL_NPZ):
    print('Ex03 Fig 3: SEG salt body input velocity model')
    sm = np.load(SALT_MODEL_NPZ)
    x_km, y_km, depth_km = sm['x'] / 1e3, sm['y'] / 1e3, sm['depth'] / 1e3
    yb, db = int(sm['ybest'][0]), int(sm['dbest'][0])
    # distinct perceptually-uniform, colourblind-safe palettes per property
    VP_CMAP, VS_CMAP = 'viridis', 'magma'
    # shared colour range per property so its cross-section and map are comparable
    vp_lo, vp_hi = float(sm['VP_xz'].min()), float(sm['VP_xz'].max())
    vs_lo, vs_hi = float(sm['VS_xz'].min()), float(sm['VS_xz'].max())

    # prepend the 200 m marine water column (top region of SEG_C3.bm: Vp=1.5 km/s,
    # Vs=0) so the cross-sections show the full model from the sea surface (depth 0).
    # The sampled model data otherwise begins at the 200 m seafloor.
    WATER_VP, WATER_VS, WATER_BASE = 1.5, 0.0, 0.2          # km/s, km/s, km
    _nw, _nx = 6, sm['VP_xz'].shape[1]
    _wd = np.linspace(0.0, WATER_BASE, _nw, endpoint=False)  # 0 .. <0.2 km
    depth_km_full = np.concatenate([_wd, depth_km])
    VP_xz_full = np.vstack([np.full((_nw, _nx), WATER_VP), sm['VP_xz']])
    VS_xz_full = np.vstack([np.full((_nw, _nx), WATER_VS), sm['VS_xz']])
    vs_lo = 0.0                                              # water carries no shear

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0), constrained_layout=True)

    def _xz(ax, depth, field, cmap, lo, hi, title, clabel):
        im = ax.pcolormesh(x_km, depth, field, cmap=cmap, vmin=lo, vmax=hi,
                           shading='auto', rasterized=True)  # keep the PDF small
        ax.axhline(depth_km[db], color='w', lw=0.8, ls=':')   # the map-slice depth
        ax.axhline(WATER_BASE, color='deepskyblue', lw=1.1)   # seafloor (water base)
        ax.annotate('water (1.5 km/s)', xy=(0.0, WATER_BASE),
                    xytext=(0.0, -0.55), ha='center', va='bottom', fontsize=8,
                    arrowprops=dict(arrowstyle='->', color='deepskyblue', lw=0.9))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('x (km)'); ax.set_ylabel('Depth (km)')
        ax.invert_yaxis(); ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, shrink=0.9, label=clabel)

    def _xy(ax, field, cmap, lo, hi, title, clabel):
        # field is [x, y] -> transpose to plot (y vs x)
        im = ax.pcolormesh(x_km, y_km, field.T, cmap=cmap, vmin=lo, vmax=hi,
                           shading='auto', rasterized=True)  # keep the PDF small
        ax.axhline(y_km[yb], color='w', lw=0.8, ls='--')      # the cross-section line
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('x (km)'); ax.set_ylabel('y (km)')
        ax.set_aspect('equal')
        fig.colorbar(im, ax=ax, shrink=0.9, label=clabel)

    # top row: P-wave speed; bottom row: S-wave speed
    _xz(axes[0, 0], depth_km_full, VP_xz_full, VP_CMAP, vp_lo, vp_hi,
        f'$V_P$ cross-section (y = {y_km[yb]:.1f} km)', '$V_P$ (km/s)')
    _xy(axes[0, 1], sm['VP_xy'], VP_CMAP, vp_lo, vp_hi,
        f'$V_P$ map (depth = {depth_km[db]:.1f} km)', '$V_P$ (km/s)')
    _xz(axes[1, 0], depth_km_full, VS_xz_full, VS_CMAP, vs_lo, vs_hi,
        f'$V_S$ cross-section (y = {y_km[yb]:.1f} km)', '$V_S$ (km/s)')
    _xy(axes[1, 1], sm['VS_xy'], VS_CMAP, vs_lo, vs_hi,
        f'$V_S$ map (depth = {depth_km[db]:.1f} km)', '$V_S$ (km/s)')

    fig.suptitle('SEG/EAGE C3 salt body -- input seismic velocity model', fontsize=12)
    save_fig(fig, os.path.join(FIG_DIR, 'ex03_salt_model_velocities'), tight=False)
else:
    print(f'Ex03 Fig 3: skipped (optional input not found: {SALT_MODEL_NPZ})')

print('Ex03 figures complete.')
