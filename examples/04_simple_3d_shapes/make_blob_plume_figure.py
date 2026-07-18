#!/usr/bin/env python3
"""
Example 04: Simple 3D Shapes - publication composite of the two wavefield
time-lapses, each preceded by its input velocity model:

  TOP    row : Cartesian low-Vs blob  - input Vs model + U_Z snapshots at
               t = 6.1 s and the 3 following outputs.
  BOTTOM row : Global plume           - input V model + U_Z snapshots at
               t = 73 s and the 3 following outputs.

All colorbars sit OUTSIDE the panels (models on the left, U_Z on the right),
one per row. Large fonts; publication layout.

INPUT
-----
This figure reads four reduced NumPy archives from a `data/` directory
(`INPUT_DIR` below), NOT the raw solver output directly:

  cartesian_blob.npz          - blob input Vs model (x, z, dvs_xz)
  plume_model.npz             - plume input V model  (lat, depth, dv_xz, ...)
  cartesian_cross_section.npz - blob U_Z(s, z, t) on the phi=0 slice
  plume_cross_section.npz     - plume U_Z(s, z, t) on the phi=0 slice

The two `*_cross_section.npz` files are extracts of the element NetCDF the two
runs write under `output/elements/`; the raw element output is far larger than
the figure needs, so it is reduced to a compact npz first. See the README
"Figures" section for how to produce these archives (run the sim -> output/,
then extract the phi=0 slice). The figure is skipped with a message if any of
the four archives are missing.

Standalone (headless): `python3 make_blob_plume_figure.py`.
"""
import os
import sys
import numpy as np
import matplotlib.tri as mtri
from matplotlib.patches import Circle, Rectangle, Ellipse

# Locate everything relative to THIS file so the script is portable.
HERE = os.path.dirname(os.path.abspath(__file__))
# Make the shared pub_style helper (staged at examples/pub_style.py) importable.
sys.path.insert(0, os.path.dirname(HERE))
from pub_style import apply_style, save_fig, COLORS  # noqa: E402  (Agg set inside)
import matplotlib.pyplot as plt  # noqa: E402

# Reduced cross-section/model archives this figure reads (see module docstring).
INPUT_DIR = os.path.join(HERE, 'data')

apply_style()
plt.rcParams.update({'axes.titlesize': 15, 'axes.labelsize': 15,
                     'xtick.labelsize': 12, 'ytick.labelsize': 12})

# geometry: blob ellipse measured from the input model (>50% Vs-reduction region)
BLOB_CX, BLOB_CZ, BLOB_WX, BLOB_WZ = 10.8, 8.6, 7.8, 5.8            # km (centre, full w/h)
SRC_DEPTH = 8.0                                                      # km, on-axis source
PLUME_HEAD_DEPTH, PLUME_HEAD_R = 1200.0, 500.0                       # km
PLUME_CONDUIT_R, PLUME_CONDUIT_DTOP, PLUME_CONDUIT_DBOT = 100.0, 1000.0, 3000.0
R_TOT = 4000.0                                                       # km, plume model radius

T_BLOB, T_PLUME = 6.1, 73.0    # first snapshot of each row; +3 following outputs
# 8-snapshot grids the standalone cross-section figures sample (so this composite
# shows the SAME time steps as the cartesian / plume cross-section figures)
GRID_BLOB = (1.0, 13.0)        # np.linspace(1, 13, 8)
GRID_PLUME = (4.0, 500.0)      # np.linspace(4, 500, 8)
CB_LAB, TI_FS = 14, 14         # colorbar-label and panel-title font sizes


def four_from_grid(times, lo, hi, t0, n=8):
    """Snapshot nearest t0 on the same linspace(lo, hi, n) grid the standalone
    cross-section figures use, plus the 3 following grid points - each mapped to
    the nearest available data sample (blob: 6.1, 7.9, 9.5, 11.3 s;
    plume: 73, 148, 217, 287 s)."""
    grid = np.linspace(lo, hi, n)
    g0 = min(int(np.argmin(np.abs(grid - t0))), n - 4)
    return [int(np.argmin(np.abs(times - grid[g]))) for g in range(g0, g0 + 4)]


def wave_panels(npz):
    z = np.load(npz)
    coords, wave, times = z['coords_sz'], z['wave_zt'], z['times']
    s_km = coords[:, 0] / 1e3
    depth_km = (coords[:, 1].max() - coords[:, 1]) / 1e3
    tri = mtri.Triangulation(s_km, depth_km)
    return tri, wave, times, s_km, depth_km


def main():
    blob_npz = os.path.join(INPUT_DIR, 'cartesian_cross_section.npz')
    plume_npz = os.path.join(INPUT_DIR, 'plume_cross_section.npz')
    bmod_npz = os.path.join(INPUT_DIR, 'cartesian_blob.npz')
    pmod_npz = os.path.join(INPUT_DIR, 'plume_model.npz')
    for f in (blob_npz, plume_npz, bmod_npz, pmod_npz):
        if not os.path.isfile(f):
            print(f'  [skip] blob/plume composite: missing {os.path.basename(f)} '
                  f'(expected under {INPUT_DIR}; see README "Figures").')
            return

    btri, bwave, btimes, _, bdep = wave_panels(blob_npz)
    ptri, pwave, ptimes, _, pdep = wave_panels(plume_npz)
    bsnaps = four_from_grid(btimes, *GRID_BLOB, T_BLOB)
    psnaps = four_from_grid(ptimes, *GRID_PLUME, T_PLUME)

    # per-row diverging colour scale for U_Z (percentile clips the near-source spike)
    bvmax = np.percentile(np.abs(bwave[:, bsnaps]), 99.5) or 1e-30
    pvmax = np.percentile(np.abs(pwave[:, psnaps]), 99.5) or 1e-30
    blev = np.linspace(-bvmax, bvmax, 41)
    plev = np.linspace(-pvmax, pvmax, 41)

    # framing per row (distance s, depth)
    BX, BD = 18.0, 17.0          # blob: km
    PX, PD = 1600.0, 3200.0      # plume: km - full head + conduit + downgoing front

    fig = plt.figure(figsize=(19.5, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 5, wspace=0.03)
    axb = [fig.add_subplot(gs[0, c]) for c in range(5)]
    axp = [fig.add_subplot(gs[1, c]) for c in range(5)]

    # ---------- TOP ROW: blob ----------
    bz = np.load(bmod_npz)
    xk = bz['x'] / 1e3
    zk = bz['z'] / 1e3
    dvs = -bz['dvs_xz'] * 100.0          # % velocity (Vp) reduction, (x, z)
    pos = xk >= 0
    im_bm = axb[0].pcolormesh(xk[pos], zk, dvs[pos].T, cmap='viridis',
                              vmin=0, vmax=100, shading='auto')
    axb[0].set_title('Input $V_P$ model', fontsize=TI_FS)
    cf_b = None
    for k, ti in enumerate(bsnaps):
        ax = axb[k + 1]
        cf_b = ax.tricontourf(btri, bwave[:, ti], levels=blev, cmap='RdBu_r', extend='both')
        cf_b.set_rasterized(True)
        ax.set_title(f't = {btimes[ti]:.1f} s', fontsize=TI_FS)
    for ax in axb:
        ax.add_patch(Ellipse((BLOB_CX, BLOB_CZ), BLOB_WX, BLOB_WZ, fill=False,
                             ec=COLORS['tertiary'], lw=1.3, ls='--', zorder=5))
        ax.plot(0.0, SRC_DEPTH, '*', ms=22, color=COLORS['secondary'],
                mec='k', mew=0.8, zorder=6)
        ax.set_aspect('equal'); ax.set_xlim(0, BX); ax.set_ylim(BD, 0)
        ax.tick_params(labelsize=12)
    axb[0].set_ylabel('Depth (km)')
    for ax in axb[1:]:
        ax.set_yticklabels([])

    # ---------- BOTTOM ROW: plume ----------
    # Reproject the input model into the SAME geometry as the wavefield panels:
    # depth measured from the top of the spherical slice, so the free surface shows
    # as the curved arc (the "circle (surface)") instead of a flat top.
    pm = np.load(pmod_npz)
    lat = pm['lat']; depth_km = pm['depth'] / 1e3
    LAT, DEP = np.meshgrid(np.radians(lat), depth_km, indexing='ij')
    r = R_TOT - DEP
    S = r * np.sin(LAT)
    ZSURF = R_TOT - r * np.cos(LAT)          # depth below top-of-slice (matches wavefield)
    red_xz = -pm['dv_xz'] * 100.0
    im_pm = axp[0].pcolormesh(S, ZSURF, red_xz, cmap='viridis', vmin=0, vmax=10,
                              shading='auto')
    axp[0].set_title('Input $V_P$ model', fontsize=TI_FS)
    cf_p = None
    for k, ti in enumerate(psnaps):
        ax = axp[k + 1]
        cf_p = ax.tricontourf(ptri, pwave[:, ti], levels=plev, cmap='RdBu_r', extend='both')
        cf_p.set_rasterized(True)
        ax.set_title(f't = {ptimes[ti]:.0f} s', fontsize=TI_FS)
    # free-surface arc (Earth's surface at radius R_TOT) - drawn on every plume panel
    s_arc = np.linspace(0.0, PX, 240)
    d_arc = R_TOT - np.sqrt(np.maximum(R_TOT ** 2 - s_arc ** 2, 0.0))
    for ax in axp:
        ax.plot(s_arc, d_arc, color=COLORS['gray'], lw=1.4, zorder=6)
        ax.add_patch(Circle((0, PLUME_HEAD_DEPTH), PLUME_HEAD_R, fill=False,
                            ec=COLORS['tertiary'], lw=1.1, ls='--', zorder=5))
        ax.add_patch(Rectangle((0, PLUME_CONDUIT_DTOP), PLUME_CONDUIT_R,
                               PLUME_CONDUIT_DBOT - PLUME_CONDUIT_DTOP, fill=False,
                               ec=COLORS['tertiary'], lw=1.1, ls='--', zorder=5))
        ax.set_aspect('equal'); ax.set_xlim(0, PX); ax.set_ylim(PD, 0)
        ax.set_xlabel('Distance (km)'); ax.tick_params(labelsize=12)
    axp[0].set_ylabel('Depth (km)')
    for ax in axp[1:]:
        ax.set_yticklabels([])

    # ---------- colorbars, all OUTSIDE the panels ----------
    # attach both colorbars to the WHOLE row (not the model alone / snapshots alone)
    # so the five panels stay one evenly-spaced group with no gap; colorbars sit at
    # the row's outer left/right edges.
    cb = fig.colorbar(im_bm, ax=axb, location='left', shrink=0.85, pad=0.01)
    cb.set_label('$V_P$ reduction (%)', fontsize=CB_LAB)
    cb_uzb = fig.colorbar(cf_b, ax=axb, location='right', shrink=0.85, pad=0.01)
    cb_uzb.set_label('$U_Z$ (m)', fontsize=CB_LAB)
    cb = fig.colorbar(im_pm, ax=axp, location='left', shrink=0.85, pad=0.01)
    cb.set_label('$V_P$ reduction (%)', fontsize=CB_LAB)
    cb_uzp = fig.colorbar(cf_p, ax=axp, location='right', shrink=0.85, pad=0.01)
    cb_uzp.set_label('$U_Z$ (m)', fontsize=CB_LAB)
    # readable U_Z ticks: a few round values shown with a scientific offset factor
    # (x10^-3 for the blob, x10^-11 for the plume) instead of long 0.000959 decimals
    import matplotlib.ticker as mticker
    for cbu, vm in ((cb_uzb, bvmax), (cb_uzp, pvmax)):
        ticks = mticker.MaxNLocator(nbins=5, symmetric=True).tick_values(-vm, vm)
        cbu.set_ticks(ticks[(ticks >= -vm) & (ticks <= vm)])
        sf = mticker.ScalarFormatter(useMathText=True)
        sf.set_powerlimits((-2, 3))
        cbu.ax.yaxis.set_major_formatter(sf)
        cbu.ax.yaxis.get_offset_text().set_fontsize(12)
        cbu.ax.tick_params(labelsize=12)

    # shape-outline / source legend (outside, top of figure)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=COLORS['tertiary'], lw=1.2, ls='--',
                      label='Velocity anomaly (blob / plume)'),
               Line2D([0], [0], color=COLORS['secondary'], marker='*', lw=0, ms=18,
                      mec='k', mew=0.8, label='Source'),
               Line2D([0], [0], color=COLORS['gray'], lw=1.4,
                      label='Free surface (plume)')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.055),
               ncol=3, frameon=False, fontsize=13)

    # publication-style (a)/(b) row headers at the top-left of each row
    axb[0].annotate('(a) Cartesian simulation - spherical low-velocity blob',
                    xy=(0.0, 1.16), xycoords='axes fraction', fontsize=18,
                    fontweight='bold', va='bottom', ha='left', annotation_clip=False)
    axp[0].annotate('(b) Spherical (global) simulation - mantle plume',
                    xy=(0.0, 1.16), xycoords='axes fraction', fontsize=18,
                    fontweight='bold', va='bottom', ha='left', annotation_clip=False)
    save_fig(fig, os.path.join(HERE, 'figures', 'ex04_blob_plume_snapshots'), tight=True)
    print('  done: figures/ex04_blob_plume_snapshots')


if __name__ == '__main__':
    main()
