#!/usr/bin/env python3
"""
Example 05: Global anisotropy — SKS shear-wave-splitting "sticks" map.

Reproduces the splitting-comparison panel of the anisotropy paper figure: at each
station a fast-axis "stick" is drawn, oriented along the measured fast polarization
azimuth phi and with length proportional to the split delay time dt. Two measurements
are overlaid per station for the 2012-07-03 New Zealand deep event recorded by the US
Transportable Array (the "US32" line):

  * RED  = real data, Long et al. (2009) SKS splitting (Long_splitting.txt)
  * BLUE = synthetic, AxiSEM3D US32 anisotropic run (20120703_phi_dt_SKS*.txt)

Both sticks are centred on the station and use the SAME length scale, so the delay
times are directly comparable (the legend doubles as a 1 s reference). A small
orthographic inset shows the source-receiver geometry.

INPUT: the small per-station splitting measurements shipped in data/splitting/. These
are post-measurement results (not raw seismograms); they ship with the example because
they are tiny. Standalone: run ``python3 make_splitting_figure.py`` from this folder.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # examples/ root holds pub_style.py

from pub_style import apply_style, save_fig, COLORS
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

apply_style()
# manual axes placement + insets: constrained_layout would fight us here
plt.rcParams['figure.constrained_layout.use'] = False

SPLIT_DIR = os.path.join(HERE, 'data', 'splitting')

# --- map region (Pacific Northwest, the US32 line) ---
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = -122.3, -115.3, 40.8, 47.2
LAT0 = 0.5 * (LAT_MIN + LAT_MAX)

# Full stick length (km) for a 1 s delay. Tuned so dt~2-3 s sticks span a legible
# fraction of the ~560 km-wide map without saturating it.
KM_PER_S = 48.0
LW = 1.5                # stick line width

REAL_FILE = 'Long_splitting.txt'                 # Long et al. (2009)
SYN_FILE = '20120703_phi_dt_SKS_vel.txt'         # AxiSEM3D US32 synthetic
SOURCE_LAT, SOURCE_LON = -40.023, 173.756        # 2012-07-03 NZ event

LAND = '#efe6cf'        # pale tan land fill (cornsilk-like)
COAST = COLORS['gray']
BORDER = '#3a3a3a'      # dark grey — state/country lines must stay clearly visible
C_REAL = COLORS['secondary']    # vermillion  (colourblind-safe red)
C_SYN = COLORS['primary']       # blue


def load_split(fname, header):
    """Return lon, lat, phi(deg, azimuth CW from N), dt(s)."""
    p = os.path.join(SPLIT_DIR, fname)
    d = np.genfromtxt(p, skip_header=header, usecols=(0, 1, 2, 3))
    d = d[~np.isnan(d).any(axis=1)]
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def stick_segments(lon, lat, phi, dt, geod):
    """Geodesic two-headed sticks centred on each station: a line from azimuth phi to
    phi+180, total length dt * KM_PER_S (geographically correct on the projection
    regardless of latitude)."""
    pts = np.column_stack([lon, lat]).astype(float)
    half_m = 0.5 * dt * KM_PER_S * 1000.0
    a = np.asarray(geod.direct(pts, phi, half_m))           # (N,3) lon,lat,az
    b = np.asarray(geod.direct(pts, phi + 180.0, half_m))
    return [[(b[i, 0], b[i, 1]), (a[i, 0], a[i, 1])] for i in range(len(lon))]


def draw_map(fig, ax, ax_rect, *, label_states=True, reference=True, inset=True):
    """Draw the full splitting-sticks map onto `ax` (a Mercator GeoAxes whose
    figure-fraction rectangle is `ax_rect`). Reusable by both the standalone figure and
    the composite; the colour legend is added by the caller (so it can live outside the
    axes)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.geodesic import Geodesic
    import matplotlib.ticker as mticker
    pc = ccrs.PlateCarree()
    geod = Geodesic()

    rlon, rlat, rphi, rdt = load_split(REAL_FILE, 0)
    slon, slat, sphi, sdt = load_split(SYN_FILE, 2)

    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=pc)

    # --- basemap ---
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=LAND, zorder=0)
    ax.add_feature(cfeature.LAKES.with_scale('50m'), facecolor='white',
                   edgecolor=COAST, linewidth=0.3, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), edgecolor=COAST,
                   linewidth=0.7, zorder=6)
    # state + country boundaries: above the sticks and thick/dark so they read clearly
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor=BORDER,
                   linewidth=1.3, zorder=6)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor=BORDER,
                   linewidth=1.6, zorder=6)

    gl = ax.gridlines(crs=pc, draw_labels=True, linewidth=0.3,
                      color='0.8', linestyle=':', zorder=1)
    gl.top_labels = gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-122, -120, -118, -116])
    gl.ylocator = mticker.FixedLocator([42, 44, 46])
    gl.xlabel_style = {'size': 13}
    gl.ylabel_style = {'size': 13}

    # --- splitting sticks (real then synthetic, both to the same scale) ---
    ax.add_collection(LineCollection(
        stick_segments(rlon, rlat, rphi, rdt, geod), transform=pc,
        colors=C_REAL, linewidths=LW, capstyle='round', zorder=4))
    ax.add_collection(LineCollection(
        stick_segments(slon, slat, sphi, sdt, geod), transform=pc,
        colors=C_SYN, linewidths=LW, capstyle='round', zorder=5))

    # --- state labels ---
    if label_states:
        for name, lo, la in [('WASHINGTON', -118.5, 46.95), ('OREGON', -120.7, 42.9),
                             ('IDAHO', -115.9, 45.4), ('NEVADA', -116.6, 41.3),
                             ('CALIFORNIA', -121.6, 41.25)]:
            ax.text(lo, la, name, transform=pc, fontsize=13, color='#33405e',
                    style='italic', ha='center', va='center', zorder=7,
                    path_effects=_halo())

    if reference:
        _add_reference(ax)
    if inset:
        _add_globe_inset(fig, ax_rect, ccrs, cfeature, pc, rlon, rlat)


def main():
    try:
        import cartopy.crs as ccrs
    except ImportError:
        print('  (cartopy missing, skipping splitting-sticks figure)')
        return

    if not os.path.isfile(os.path.join(SPLIT_DIR, REAL_FILE)):
        print(f'  [skip] splitting data not found in {SPLIT_DIR}')
        return

    proj = ccrs.Mercator(central_longitude=0.5 * (LON_MIN + LON_MAX),
                         min_latitude=LAT_MIN, max_latitude=LAT_MAX)

    fig = plt.figure(figsize=(7.4, 8.6))
    ax_rect = [0.085, 0.165, 0.88, 0.755]
    ax = fig.add_axes(ax_rect, projection=proj)
    draw_map(fig, ax, ax_rect, reference=False)
    _add_legend(ax, 7.4, ax_rect[2])

    fig.suptitle('SKS shear-wave splitting — Pacific Northwest (US32 line)',
                 fontsize=16, y=0.95)
    save_fig(fig, os.path.join(HERE, 'figures', 'ex05_splitting_sticks'), tight=False)
    print('  done: ex05_splitting_sticks')


def _halo():
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=2.2, foreground='white')]


def _add_reference(ax):
    """A 1 s reference stick (both colours) + north arrow in a white card, upper-left
    corner. The bar length is the to-scale on-map length of a 1 s delay at the map's
    mid-latitude."""
    from matplotlib.patches import FancyBboxPatch
    frac = (KM_PER_S / (111.32 * np.cos(np.radians(LAT0)))) / (LON_MAX - LON_MIN)
    card = FancyBboxPatch((0.018, 0.852), 0.33, 0.135, transform=ax.transAxes,
                          boxstyle='round,pad=0.006,rounding_size=0.012',
                          facecolor='white', edgecolor='0.6', linewidth=0.6,
                          zorder=8)
    ax.add_patch(card)
    xc = 0.125
    for y, c in [(0.935, COLORS['secondary']), (0.890, COLORS['primary'])]:
        ax.plot([xc - frac / 2, xc + frac / 2], [y, y], transform=ax.transAxes,
                color=c, lw=LW + 0.6, solid_capstyle='round', zorder=9)
    ax.text(xc, 0.968, '1 s splitting', transform=ax.transAxes, fontsize=12,
            ha='center', va='center', zorder=9)
    ax.annotate('N', xy=(0.295, 0.965), xytext=(0.295, 0.880),
                xycoords='axes fraction', ha='center', va='center',
                fontsize=13, fontweight='bold',
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.2),
                zorder=10)


def legend_handles():
    """Shared colour-key handles (real vs synthetic) for an external legend."""
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], color=COLORS['secondary'], lw=LW + 0.6,
                   solid_capstyle='round', label='Real data — Long et al. (2009)'),
            Line2D([0], [0], color=COLORS['primary'], lw=LW + 0.6,
                   solid_capstyle='round', label='Synthetic — AxiSEM3D (US32)')]


def one_second_handlelength(fig_width_in, ax_width_frac, fontsize):
    """Legend `handlelength` (in font-size units) so each sample line equals the on-map
    length of a 1 s splitting stick (E-W at the map mid-latitude). This lets the
    colour-key legend double as the splitting scale bar."""
    frac = (KM_PER_S / (111.32 * np.cos(np.radians(LAT0)))) / (LON_MAX - LON_MIN)
    return frac * ax_width_frac * fig_width_in * 72.0 / fontsize


def _add_legend(ax, fig_width_in, ax_width_frac, fontsize=12.5):
    hl = one_second_handlelength(fig_width_in, ax_width_frac, fontsize)
    leg = ax.legend(handles=legend_handles(), loc='lower center',
                    bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False,
                    fontsize=fontsize, handlelength=hl, columnspacing=1.0,
                    handletextpad=0.5, labelspacing=0.3, borderpad=0.2,
                    title='Stick length = 1 s splitting')
    leg.get_title().set_fontsize(fontsize)


def _add_globe_inset(fig, ax_rect, ccrs, cfeature, pc, rlon, rlat):
    """Small orthographic globe showing the NZ source -> US array geometry, placed in
    the lower-left corner of the map axes (`ax_rect` = the map's figure-fraction
    rectangle, so the inset travels with the panel)."""
    def _unit(la, lo):
        la, lo = np.radians(la), np.radians(lo)
        return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    clat0 = float(np.mean(rlat)); clon0 = float(np.mean(rlon))
    mid = _unit(SOURCE_LAT, SOURCE_LON) + _unit(clat0, clon0)
    clon = np.degrees(np.arctan2(mid[1], mid[0]))
    clat = np.degrees(np.arctan2(mid[2], np.hypot(mid[0], mid[1])))

    from matplotlib.patches import FancyBboxPatch
    x0, y0, w, h = ax_rect
    card_rect = [x0 + 0.074 * w, y0 + 0.031 * h, 0.267 * w, 0.283 * h]
    bgax = fig.add_axes(card_rect, zorder=20)
    bgax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, transform=bgax.transAxes,
                                  boxstyle='round,pad=0.0,rounding_size=0.04',
                                  facecolor='white', edgecolor='0.6', linewidth=0.6))
    bgax.set_axis_off()
    gax = fig.add_axes([x0 + 0.0875 * w, y0 + 0.0415 * h, 0.2386 * w, 0.2579 * h],
                       projection=ccrs.Orthographic(clon, clat), zorder=21)
    gax.set_facecolor('white')
    gax.set_global()
    gax.add_feature(cfeature.LAND.with_scale('110m'), facecolor=LAND, zorder=0)
    gax.add_feature(cfeature.COASTLINE.with_scale('110m'), edgecolor=COAST,
                    linewidth=0.3, zorder=1)
    gax.gridlines(linewidth=0.25, color='0.7', linestyle=':')
    gax.plot([SOURCE_LON, clon0], [SOURCE_LAT, clat0], transform=ccrs.Geodetic(),
             color='0.35', lw=0.8, zorder=2)
    gax.scatter([clon0], [clat0], transform=pc, marker='s', s=22,
                c=COLORS['primary'], edgecolors='k', linewidths=0.3, zorder=4)
    gax.scatter([SOURCE_LON], [SOURCE_LAT], transform=pc, marker='*', s=150,
                c=COLORS['secondary'], edgecolors='k', linewidths=0.4, zorder=5)
    gax.set_title('source → array', fontsize=11, pad=3)
    gax.spines['geo'].set_linewidth(0.5)


if __name__ == '__main__':
    main()
