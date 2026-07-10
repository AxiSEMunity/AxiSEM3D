#!/usr/bin/env python3
"""
Adrian kernel 1-D example — publication figures (re-styled from the axikernels
HDF5 output, so they match the rest of the example gallery).
  Fig 1: P-wave (vp) cross-correlation traveltime sensitivity kernel — the
         banana-doughnut on the source-receiver great-circle slice.
  Fig 2: Moho interface sensitivity kernel (K_SS) map.

Reads kernel_output/vp_kernel.h5 (pandas: 'df' + 'metadata') and
kernel_output/moho_kd.h5 (h5py: lat_deg, lon_deg, Kd, K_dn, K_dv).
Self-contained: imports only pub_style. Run with the axikernels conda python.
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# shared pub_style.py lives at the examples/ root, two levels up
# (examples/07_kernels/kernel_1D/ -> examples/07_kernels/ -> examples/)
EXAMPLES_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, EXAMPLES_ROOT)
from pub_style import apply_style, save_fig, COLORS
import matplotlib.pyplot as plt

apply_style()
# INPUT: directory holding the kernel HDF5 that compute_kernels.py writes
# (vp_kernel.h5, moho_kd.h5). Defaults to ./kernel_output/ produced by running
# this example end to end. Override via the KERNEL_OUTPUT environment variable.
INPUT = os.environ.get('KERNEL_OUTPUT', os.path.join(HERE, 'kernel_output'))
KO = INPUT
FIG_DIR = os.path.join(HERE, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

R_SURF = 6371.0      # km
R_CMB = 3480.0
TITLE = '1-D PREM: $v_P$ traveltime sensitivity kernel (P400P)'


def plot_vp_kernel():
    import pandas as pd
    f = os.path.join(KO, 'vp_kernel.h5')
    if not os.path.exists(f):
        print('  vp_kernel.h5 missing, skipping'); return
    df = pd.read_hdf(f, 'df')
    md = pd.read_hdf(f, 'metadata')
    # slice is in the lat=0 plane -> Cartesian in-plane coords from (radius, lon)
    r = df['radius'].values / 1e3
    lon = df['lon'].values
    x = r * np.cos(lon)
    y = r * np.sin(lon)
    data = df['data'].values
    # The kernel spans a huge dynamic range — a near-receiver spike is ~100x the
    # banana-doughnut body, so the axikernels 0.1*max scale washes the structure
    # out to white. Clip at the 99th percentile of nonzero |K| so the doughnut is
    # visible; the spike saturates (which is fine).
    absd = np.abs(data[np.isfinite(data) & (data != 0)])
    vmax = float(np.percentile(absd, 99)) if absd.size else (np.nanmax(np.abs(data)) or 1.0)
    if not np.isfinite(vmax) or vmax == 0:
        vmax = np.nanmax(np.abs(data)) or 1.0

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    tcf = ax.tricontourf(x, y, data, levels=np.linspace(-vmax, vmax, 41),
                         cmap='RdBu_r', extend='both')
    cb = fig.colorbar(tcf, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label('$K_{v_P}$ (s m$^{-3}$)', fontsize=9)
    # reference interfaces
    th = np.linspace(-np.pi, np.pi, 400)
    for rr, lab in [(R_SURF, None), (R_CMB, None)]:
        ax.plot(rr * np.cos(th), rr * np.sin(th), color=COLORS['gray'], lw=0.6)
    # source (point1) and receiver (point2) from metadata
    p1 = np.array(md['point1']) / 1e3
    p2 = np.array(md['point2']) / 1e3
    ax.plot(p1[0], p1[1], '*', ms=16, color='k', mec='w', mew=0.5, zorder=10, label='source')
    ax.plot(p2[0], p2[1], '^', ms=10, color=COLORS['tertiary'], mec='k', mew=0.4,
            zorder=10, label='receiver')
    ax.set_aspect('equal')
    ax.set_xlabel('x (km)'); ax.set_ylabel('y (km)')
    ax.set_title(TITLE, fontsize=11)
    ax.legend(loc='lower left', fontsize=8, frameon=True)
    save_fig(fig, os.path.join(FIG_DIR, 'kernel1D_vp_slice'))


def plot_moho_kernel(topo_nc=None, title='1-D PREM: Moho interface kernel $K_{SS}$'):
    import h5py
    f = os.path.join(KO, 'moho_kd.h5')
    if not os.path.exists(f):
        print('  moho_kd.h5 missing, skipping'); return
    with h5py.File(f, 'r') as h:
        lat = h['lat_deg'][:]; lon = h['lon_deg'][:]; Kd = h['Kd'][:]
    finite = np.abs(Kd[np.isfinite(Kd) & (Kd != 0)])
    vmax = np.percentile(finite, 99) if finite.size else (np.nanmax(np.abs(Kd)) or 1.0)
    fig, ax = plt.subplots(figsize=(7, 4))
    pm = ax.pcolormesh(lon, lat, Kd, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                       shading='auto', rasterized=True)
    cb = fig.colorbar(pm, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label('$K_{SS}$', fontsize=9)
    if topo_nc and os.path.exists(topo_nc):
        try:
            import netCDF4
            d = netCDF4.Dataset(topo_nc)
            tlat = d.variables['latitude'][:]; tlon = d.variables['longitude'][:]
            und = d.variables['undulation_MOHO'][:]
            TL, TA = np.meshgrid(tlon, tlat)
            ax.contour(TL, TA, und, colors='k', linewidths=0.4, alpha=0.6)
            d.close()
        except Exception as e:
            print('  topo overlay skipped:', e)
    ax.plot(0, 0, '*', ms=16, color='k', mec='w', mew=0.5, zorder=10)
    ax.plot(40, 0, '^', ms=10, color=COLORS['tertiary'], mec='k', mew=0.4, zorder=10)
    ax.set_xlim(float(np.min(lon)), float(np.max(lon)))
    ax.set_ylim(float(np.min(lat)), float(np.max(lat)))
    ax.set_xlabel('Longitude (deg)'); ax.set_ylabel('Latitude (deg)')
    ax.set_title(title, fontsize=11)
    save_fig(fig, os.path.join(FIG_DIR, 'kernel1D_moho'))


if __name__ == '__main__':
    print('Kernel 1D Fig 1: vp slice'); plot_vp_kernel()
    print('Kernel 1D Fig 2: Moho kernel'); plot_moho_kernel()
    print('Kernel 1D figures complete.')
