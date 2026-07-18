"""
Shared publication-quality plotting style for AxiSEM3D examples.
Import this module before creating any figures.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

STYLE = {
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'lines.linewidth': 1.0,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'mathtext.fontset': 'cm',
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.constrained_layout.use': True,
}

# Colorblind-safe palette (Okabe-Ito). Avoids red/green confusion.
COLORS = {
    'primary': '#0072B2',     # blue
    'secondary': '#D55E00',   # vermillion
    'tertiary': '#009E73',    # bluish green
    'quaternary': '#E69F00',  # orange
    'quinary': '#CC79A7',     # reddish purple
    'accent': '#009E73',
    'R': '#D55E00',           # vermillion
    'T': '#009E73',           # bluish green
    'Z': '#0072B2',           # blue
    'mono': '#000000',        # single-category traces (record sections)
    'gray': '#636363',
    'light_gray': '#bdbdbd',
}

CHANNEL_COLORS = {'R': COLORS['R'], 'T': COLORS['T'], 'Z': COLORS['Z']}
CHANNEL_LABELS = {'R': 'Radial', 'T': 'Transverse', 'Z': 'Vertical'}


def apply_style():
    plt.rcParams.update(STYLE)


def epicentral_distance_deg(src_lat, src_lon, st_lat, st_lon):
    """Great-circle angular distance source->station, in degrees."""
    import numpy as np
    la1, lo1, la2, lo2 = map(np.radians, (src_lat, src_lon, st_lat, st_lon))
    cosd = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    return float(np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0))))


def save_fig(fig, path_stem, tight=True):
    """Save as both PDF and PNG."""
    p = Path(path_stem)
    p.parent.mkdir(parents=True, exist_ok=True)
    kwargs = {'bbox_inches': 'tight'} if tight else {}
    fig.savefig(str(p.with_suffix('.pdf')), **kwargs)
    fig.savefig(str(p.with_suffix('.png')), **kwargs)
    print(f'  Saved {p.with_suffix(".pdf")} and {p.with_suffix(".png")}')
    plt.close(fig)


def bandpass(data, dt, period_min, period_max=None, corners=4, axis=-1):
    """Zero-phase Butterworth bandpass for seismograms.

    `period_min` is the SHORT-period (high-frequency) corner: set it near the
    mesh's resolved/dominant period (e.g. 5 s for a 5 s run) — the solver has no
    reliable energy below it. `period_max` is the LONG-period (low-frequency)
    corner and can be generous (e.g. a large fraction of the record length);
    pass None for a high-pass-only (low-cut) filter. Filtering is applied along
    `axis` with filtfilt (zero phase, no time shift)."""
    import numpy as np
    from scipy.signal import butter, filtfilt
    fnyq = 0.5 / dt
    f_hi = min((1.0 / period_min) / fnyq, 0.99)        # short period -> high freq
    if period_max is None or not np.isfinite(period_max):
        b, a = butter(corners, f_hi, btype='lowpass')
    else:
        f_lo = max((1.0 / period_max) / fnyq, 1e-6)    # long period -> low freq
        b, a = butter(corners, [f_lo, f_hi], btype='bandpass')
    return filtfilt(b, a, data, axis=axis)


def source_receiver_map(src_latlon, recv_latlon, path_stem, *, title='',
                        coastlines=True, src_label='source', recv_label=None,
                        body='Earth'):
    """Source-receiver configuration map (ex00 style: white background, coarse
    grey coastlines only where applicable, orange star = source, blue triangles =
    receivers) with great-circle source-receiver paths drawn.

    src_latlon = (lat, lon) of the source; recv_latlon = list of (lat, lon).
    coastlines=True draws Earth coastlines; set False for other bodies (Mars,
    Moon) where the orthographic globe outline alone frames the geometry.
    The orthographic is centred on the source<->receiver-cluster midpoint so a
    single source and a tight receiver array are both well inside the disk."""
    import numpy as np
    recv = [r for r in recv_latlon if r is not None]
    if not recv:
        return
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        print('  (cartopy missing, skipping source-receiver map)')
        return
    rlat = np.array([r[0] for r in recv]); rlon = np.array([r[1] for r in recv])
    slat, slon = src_latlon
    if recv_label is None:
        recv_label = f'receivers (n={len(recv)})'

    def _unit(la, lo):
        la, lo = np.radians(la), np.radians(lo)
        return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    # centre on the great-circle MIDPOINT between the source and the FARTHEST
    # receiver, so the whole span fits inside the visible hemisphere (a long
    # meridian profile, source->135 deg, would otherwise push the far stations
    # off/onto the limb of an orthographic centred on the receiver centroid).
    su = _unit(slat, slon)
    ru = _unit(rlat, rlon)
    ifar = int(np.argmin((su[:, None] * ru).sum(axis=0)))
    mid = su + _unit(rlat[ifar], rlon[ifar])
    clon = np.degrees(np.arctan2(mid[1], mid[0]))
    clat = np.degrees(np.arctan2(mid[2], np.hypot(mid[0], mid[1])))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(clon, clat))
    ax.set_global()
    if coastlines:
        ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4,
                       color=COLORS['gray'])
    for la_r, lo_r in zip(rlat, rlon):
        ax.plot([slon, lo_r], [slat, la_r], transform=ccrs.Geodetic(),
                color=COLORS['primary'], lw=0.3, alpha=0.35, zorder=3)
    ax.scatter(rlon, rlat, transform=ccrs.PlateCarree(), s=24, c=COLORS['primary'],
               marker='^', edgecolors='k', linewidths=0.3, zorder=10, label=recv_label)
    ax.scatter([slon], [slat], transform=ccrs.PlateCarree(), s=240, c=COLORS['secondary'],
               marker='*', edgecolors='k', linewidths=0.5, zorder=20, label=src_label)
    ax.legend(loc='lower left', frameon=True, framealpha=0.9, fontsize=8)
    ax.set_title(title, fontsize=11)
    save_fig(fig, path_stem, tight=False)


def sequential_cmap(n):
    """Return n evenly-spaced colors from a perceptually uniform colormap."""
    cmap = plt.cm.viridis
    return [cmap(i / max(1, n - 1)) for i in range(n)]


def active_window_indices(data, n, thresh_frac=0.04):
    """Pick n evenly-spaced time indices over the window where the array
    amplitude exceeds thresh_frac of its peak (so snapshots land on the
    actual wave-propagation window, not flat pre/post-signal time)."""
    import numpy as np
    amp = np.max(np.abs(data), axis=1)
    peak = amp.max() if amp.size else 0.0
    active = np.where(amp > thresh_frac * peak)[0] if peak > 0 else np.array([])
    if active.size >= n:
        i0, i1 = int(active[0]), int(active[-1])
    else:
        i0, i1 = 0, len(data) - 1
    return np.linspace(i0, i1, n).astype(int)


def wavefield_timelapse_map(lons, lats, data, times, path_stem, *,
                            nrows=2, ncols=4, cmap='RdBu_r', clabel='',
                            title='', aspect=1.3, vmax=None, pct=99.5,
                            t_range=None):
    """Static time-lapse: an nrows x ncols grid of map snapshots of a scalar
    wavefield (one panel per time), shared colour scale + one colorbar.
    Replaces a movie for the manuscript main text.
      data : (ntime, nstation)   lons,lats : (nstation,)   times : (ntime,)
      t_range=(tmin,tmax): pick the snapshots evenly across this time window
        (seconds); if None, use the auto amplitude-active window.
    """
    import numpy as np
    n = nrows * ncols
    if t_range is not None:
        target = np.linspace(t_range[0], t_range[1], n)
        snap = np.array([int(np.argmin(np.abs(times - tt))) for tt in target])
    else:
        snap = active_window_indices(data, n)
    if vmax is None:
        vmax = np.percentile(np.abs(data[snap]), pct) or 1e-30
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.7 * ncols, 2.7 * nrows),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    sc = None
    for k, ti in enumerate(snap):
        ax = axes[k]
        sc = ax.scatter(lons, lats, s=2, c=data[ti], vmin=-vmax, vmax=vmax,
                        cmap=cmap, rasterized=True)
        ax.set_title(f't = {times[ti]:.0f} s', fontsize=9)
        ax.set_aspect(aspect)
        ax.tick_params(labelsize=7)
        if k % ncols:
            ax.set_yticklabels([])
        if k < n - ncols:
            ax.set_xticklabels([])
    for k in range(len(snap), len(axes)):
        axes[k].set_visible(False)
    cbar = fig.colorbar(sc, ax=axes.tolist(), shrink=0.6, pad=0.02)
    cbar.set_label(clabel, fontsize=9)
    if title:
        fig.suptitle(title, fontsize=12)
    save_fig(fig, path_stem, tight=False)
