#!/usr/bin/env python3
"""
Example 05: Global full (3-D) anisotropy — SKS / SKKS record sections.

Stations on the meridian cases lie along the source meridian spanning epicentral
distance Delta = 85-135 deg (the SKS window). The source is at the north pole, so
AxiSEM3D's RTZ station output gives the radial and transverse displacement directly
(channels U in RTZ = R, T, Z).

For each anisotropy scenario a 2-panel record section (Radial | Transverse) shows a
window around the predicted SKS arrival (aligned on SKS, PREM/taup). Radial and
transverse are normalized JOINTLY per station (same factor = joint max|R|,|T| in the
window) so the transverse energy produced by anisotropy is shown to scale relative to
the radial SKS pulse. A purely isotropic Earth would have zero transverse SKS, so the
transverse energy here is the splitting signature.

INPUT
-----
Each scenario is read from AxiSEM3D ASCII_STATION output: a ``data_time.ascii`` time
column plus one ``<network>.<station>.ascii`` file per receiver, each holding the U
displacement vector. By default this script looks for those files inside each case's
``output/stations/JW_stations/`` directory (i.e. the output you get by running the
case here), and falls back to a flat ``data/<subdir>/`` extraction if present. Set the
INPUT_ROOT variable below, or the EX05_INPUT environment variable, to point elsewhere.

IMPORTANT — resolution and coordinate frame
  The publication figures were produced from a fine (5 s) run whose station output was
  written in the RTZ frame with channels [U]. The coarse 50 s cases shipped with this
  example are a fast functional demonstration; the meridian (PREM / deep-mantle) cases
  ship with ENZ output, so their radial/transverse panels are only meaningful when re-run
  with ``coordinate_frame: RTZ``. The US32 cases already ship with RTZ output. See the
  README "Reproduce -> Figures" section for the exact steps to regenerate the publication
  figures. Scenarios with no readable RTZ output are skipped automatically.
"""
import os
import sys
import glob

import numpy as np

# --- locate inputs and the shared style helper RELATIVE to this file -------------
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # examples/ root holds pub_style.py

from pub_style import (apply_style, save_fig, COLORS, CHANNEL_LABELS,
                       epicentral_distance_deg, bandpass)
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel

apply_style()

# Root under which the seismograms are found. Defaults to this example directory; each
# scenario is resolved to its own output/ (or data/) sub-path below. Override with the
# EX05_INPUT environment variable if you keep the runs elsewhere.
INPUT_ROOT = os.environ.get('EX05_INPUT', HERE)
DATA_DIR = os.path.join(HERE, 'data')              # small shipped station/splitting metadata

MODEL = TauPyModel('prem')
WIN_PRE, WIN_POST = 15.0, 45.0
US32_EVENT = (-40.023, 173.756)                    # 2012-07-03 New Zealand event lat/lon

# Scenarios: (label, case-relative path, distance-geometry, source depth km, display_shift_s).
# The first three share the idealized north-pole source + meridian station line
# (Delta = 85-135 deg). The US32 entries are the real 2012-07-03 New Zealand deep
# earthquake recorded by the US Transportable Array, so their distances are
# great-circle from that event.
#
# display_shift_s (last field): seconds to move the plotted traces to the RIGHT, so the
# SKS *onset* (the START of the waveform), not the broadband peak, sits on the predicted
# SKS line. Only the US32 sections use it; the prediction lines themselves are not moved.
SCENARIOS = [
    ('Anisotropic PREM mesh',         'PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh',        'meridian', 0.0,   0.0),
    ('Isotropic mesh + ani. overlay', 'PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani', 'meridian', 0.0,  0.0),
    ('Lowermost-mantle anisotropy',   'deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani',     'meridian', 0.0,   0.0),
    ('US32 olivine (2012 NZ event)',  '2012-07-03_paper_example_50s/sim_US32_olivineE',                   'us32',     229.8, 5.0),
    ('US32 olivine (constant Nr)',    '2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU',            'us32',     229.8, 5.0),
]


def resolve_station_dir(case_rel):
    """Return the directory that holds the <station>.ascii files for a scenario, or None.

    Tries, in order: the case's raw run output (INPUT_ROOT/<case_rel>/output/stations/*),
    then a flat extraction layout (INPUT_ROOT/data/<basename>/ or INPUT_ROOT/<basename>/),
    so the script reads whatever this example actually produced without editing paths."""
    base = os.path.basename(case_rel)
    candidates = [
        os.path.join(INPUT_ROOT, case_rel, 'output', 'stations', 'JW_stations'),
        os.path.join(INPUT_ROOT, case_rel, 'output', 'stations'),
        os.path.join(INPUT_ROOT, 'data', base),
        os.path.join(INPUT_ROOT, base),
    ]
    for c in candidates:
        # a station group dir under output/stations/ (any single group name)
        if c.endswith('stations') and os.path.isdir(c):
            for grp in sorted(glob.glob(os.path.join(c, '*'))):
                if os.path.isfile(os.path.join(grp, 'data_time.ascii')):
                    return grp
        if os.path.isfile(os.path.join(c, 'data_time.ascii')):
            return c
    return None


def load_meridian_distances():
    """SK.D## meridian line: epicentral distance = 90 - lat (source at the pole)."""
    f = os.path.join(DATA_DIR, 'sks_stations.txt')
    name2dist = {}
    if os.path.isfile(f):
        for line in open(f):
            c = line.split('#')[0].split()
            if len(c) >= 4:
                nm, net, lat = c[0], c[1], float(c[2])
                name2dist[nm] = 90.0 - lat
                name2dist[f'{net}.{nm}'] = 90.0 - lat
    return name2dist


def load_us32_distances():
    """TA stations (us32_stations.txt): great-circle distance from the NZ event."""
    f = os.path.join(DATA_DIR, 'us32_stations.txt')
    name2dist = {}
    if os.path.isfile(f):
        for line in open(f):
            c = line.split()
            if len(c) < 4 or line.startswith('#'):
                continue
            nm, net, lat, lon = c[0], c[1], float(c[2]), float(c[3])
            d = epicentral_distance_deg(*US32_EVENT, lat, lon)
            name2dist[f'{net}.{nm}'] = d
            name2dist[nm] = d
    return name2dist


def load_meridian_coords():
    """(lat, lon) of the SK.D## meridian receivers, + the north-pole source."""
    f = os.path.join(DATA_DIR, 'sks_stations.txt')
    recv = []
    if os.path.isfile(f):
        for line in open(f):
            c = line.split('#')[0].split()
            if len(c) >= 4:
                recv.append((float(c[2]), float(c[3])))
    return (89.99999, 0.0), recv          # source at the north pole


def load_us32_coords():
    """(lat, lon) of the TA receivers, + the 2012 NZ event source."""
    f = os.path.join(DATA_DIR, 'us32_stations.txt')
    recv = []
    if os.path.isfile(f):
        for line in open(f):
            c = line.split()
            if len(c) >= 4 and not line.startswith('#'):
                recv.append((float(c[2]), float(c[3])))
    return US32_EVENT, recv


def plot_geometry(label, stem, src_latlon, recv_latlon):
    """Source-receiver configuration map (white background, coarse coastlines only,
    orange star = source, blue triangles = receivers) with the great-circle paths drawn."""
    if not recv_latlon:
        return
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
    except ImportError:
        print('  (cartopy missing, skipping geometry map)')
        return
    rlat = np.array([r[0] for r in recv_latlon])
    rlon = np.array([r[1] for r in recv_latlon])
    slat, slon = src_latlon

    def _unit(la, lo):
        la, lo = np.radians(la), np.radians(lo)
        return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    # centre on the great-circle midpoint of the source and the FARTHEST receiver so
    # the full span fits in the visible hemisphere instead of pushing far stations onto
    # the limb (handles date-line wrap, e.g. NZ -> US).
    su = _unit(slat, slon); ru = _unit(rlat, rlon)
    ifar = int(np.argmin((su[:, None] * ru).sum(axis=0)))
    mid = su + _unit(rlat[ifar], rlon[ifar])
    clon = np.degrees(np.arctan2(mid[1], mid[0]))
    clat = np.degrees(np.arctan2(mid[2], np.hypot(mid[0], mid[1])))

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(clon, clat))
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4, color=COLORS['gray'])
    for la_r, lo_r in zip(rlat, rlon):           # great-circle paths
        ax.plot([slon, lo_r], [slat, la_r], transform=ccrs.Geodetic(),
                color=COLORS['primary'], lw=0.3, alpha=0.35, zorder=3)
    ax.scatter(rlon, rlat, transform=ccrs.PlateCarree(), s=24, c=COLORS['primary'],
               marker='^', edgecolors='k', linewidths=0.3, zorder=10,
               label=f'receivers (n={len(rlat)})')
    ax.scatter([slon], [slat], transform=ccrs.PlateCarree(), s=240, c=COLORS['secondary'],
               marker='*', edgecolors='k', linewidths=0.5, zorder=20, label='source')
    ax.legend(loc='lower left', frameon=True, framealpha=0.9, fontsize=8)
    ax.set_title(f'Source-receiver geometry — {label}', fontsize=11)
    save_fig(fig, os.path.join(HERE, 'figures', f'ex05_geometry_{stem}'))


def phase_time(delta_deg, src_depth_km, phase='SKS'):
    arr = MODEL.get_travel_times(source_depth_in_km=src_depth_km,
                                 distance_in_degree=delta_deg, phase_list=[phase])
    return arr[0].time if arr else None


def measure_sks_offset(time, recs, src_depth_km, half_win=12.0):
    """Median time (s) of the observed SKS pulse relative to the taup-predicted SKS.

    For each trace take the joint envelope max(|R|, |T|) inside +/- half_win s of the
    predicted SKS and record where it peaks; return the median over stations. This is the
    empirical STF/finite-source + taup-vs-mesh alignment offset to add to t_SKS."""
    offs = []
    for dist, nm, data in recs:
        t_sks = phase_time(dist, src_depth_km, 'SKS')
        if t_sks is None:
            continue
        w = (time >= t_sks - half_win) & (time <= t_sks + half_win)
        if not w.any():
            continue
        tw = time[w] - t_sks
        env = np.maximum(np.abs(data[:, 0][w]), np.abs(data[:, 1][w]))
        if env.max() <= 0:
            continue
        offs.append(tw[int(np.argmax(env))])
    return float(np.median(offs)) if offs else 0.0


def load_scenario(station_dir, name2dist):
    if station_dir is None or not os.path.isdir(station_dir):
        return None
    tfile = os.path.join(station_dir, 'data_time.ascii')
    if not os.path.isfile(tfile):
        return None
    time = np.loadtxt(tfile)
    recs = []
    for f in sorted(glob.glob(os.path.join(station_dir, '*.ascii'))):
        nm = os.path.basename(f).replace('.ascii', '')
        if nm == 'data_time':
            continue
        dist = name2dist.get(nm)
        if dist is None:
            continue
        data = np.loadtxt(f)
        if data.ndim == 1:
            data = data[:, None]
        # Use only the first three columns = the U displacement vector (R, T, Z when the
        # case was run with coordinate_frame: RTZ). Extra channels (G, E, S, R strain/curl)
        # are ignored. We need at least two horizontal components for the R|T section.
        if data.shape[1] < 3:
            continue
        recs.append((dist, nm, data[:, :3]))
    recs.sort(key=lambda r: r[0])
    return time, recs


def plot_scenario(label, stem, station_dir, name2dist, src_depth_km=0.0, right_shift_s=0.0,
                  large_fonts=False, filt_pmin=None):
    FS = dict(title=15, panel=14, axlabel=15, tick=13, legend=11) if large_fonts \
        else dict(title=12, panel=11, axlabel=11, tick=9, legend=7.5)
    loaded = load_scenario(station_dir, name2dist)
    if loaded is None:
        print(f'  [skip] {label}: no data ({station_dir})')
        return False
    time, recs = loaded
    if not recs:
        print(f'  [skip] {label}: no matched stations')
        return False

    # Low-pass the idealized meridian runs to the SKS band (keep periods > filt_pmin s) to
    # remove unresolved high-frequency ringing; the real-event US32 runs stay unfiltered.
    if filt_pmin and len(time) > 20:
        dt = float(np.median(np.diff(time)))
        for _, _, data in recs:
            data[:, 0] = bandpass(data[:, 0], dt, filt_pmin)
            data[:, 1] = bandpass(data[:, 1], dt, filt_pmin)

    if len(recs) > 32:
        idx = np.linspace(0, len(recs) - 1, 32).astype(int)
        recs = [recs[i] for i in idx]
    span = max(recs[-1][0] - recs[0][0], 1.0)
    gain = span / 22.0

    win_post = WIN_POST
    ts_far = phase_time(recs[-1][0], src_depth_km, 'SKS')
    tk_far = phase_time(recs[-1][0], src_depth_km, 'SKKS')
    if ts_far is not None and tk_far is not None:
        win_post = float(min(max(WIN_POST, (tk_far - ts_far) + 12.0), 165.0))

    sks_offset = measure_sks_offset(time, recs, src_depth_km)
    print(f'    SKS alignment offset (observed - taup) = {sks_offset:+.2f} s')

    fig, (axR, axT) = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    dists = []
    skks_pts = []
    for dist, nm, data in recs:
        t_sks = phase_time(dist, src_depth_km, 'SKS')
        if t_sks is None:
            continue
        t_align = t_sks + sks_offset
        R = data[:, 0]
        T = data[:, 1]
        in_win = (time >= t_align - WIN_PRE - right_shift_s) & (time <= t_align + win_post)
        if not in_win.any():
            continue
        tw = time[in_win] - t_align + right_shift_s
        Rw, Tw = R[in_win], T[in_win]
        norm = max(np.max(np.abs(Rw)), np.max(np.abs(Tw)), 1e-30)
        axR.plot(tw, Rw / norm * gain + dist, color=COLORS['mono'], lw=0.6)
        axT.plot(tw, Tw / norm * gain + dist, color=COLORS['secondary'], lw=0.6)
        dists.append(dist)
        t_skks = phase_time(dist, src_depth_km, 'SKKS')
        if t_skks is not None and -WIN_PRE <= (t_skks - t_sks) <= win_post:
            skks_pts.append((t_skks - t_sks, dist))

    if not dists:
        print(f'  [skip] {label}: no traces inside the SKS window')
        plt.close(fig)
        return False

    for ax, ttl in ((axR, 'Radial'), (axT, 'Transverse')):
        ax.axvline(0.0, color=COLORS['gray'], lw=0.8, ls='--', zorder=1)
        if skks_pts:
            sp = np.array(sorted(skks_pts, key=lambda p: p[1]))
            ax.plot(sp[:, 0], sp[:, 1], color=COLORS['quaternary'], lw=1.1,
                    ls='--', zorder=2)
        ax.set_xlabel('Time relative to predicted SKS (s)', fontsize=FS['axlabel'])
        ax.set_title(ttl, fontsize=FS['panel'])
        ax.tick_params(labelsize=FS['tick'])
        ax.set_xlim(-WIN_PRE, win_post)
    axR.plot([], [], color=COLORS['gray'], lw=0.8, ls='--', label='SKS (predicted)')
    axR.plot([], [], color=COLORS['quaternary'], lw=1.1, ls='--', label='SKKS (predicted)')
    if right_shift_s:
        axR.plot([], [], color='none',
                 label=f'synthetics shifted +{right_shift_s:.0f} s\n(SKS onset on line)')
    axR.legend(loc='upper left', fontsize=FS['legend'], framealpha=0.9)
    axR.set_ylabel('Epicentral distance ($^\\circ$)', fontsize=FS['axlabel'])
    axR.set_ylim(min(dists) - gain * 1.5, max(dists) + gain * 1.5)
    band = f'low-pass {filt_pmin:.0f} s' if filt_pmin else 'broadband'
    fig.suptitle(f'SKS / SKKS record section ({band}) — {label}', fontsize=FS['title'])
    save_fig(fig, os.path.join(HERE, 'figures', f'ex05_SKS_{stem}'))
    return True


def main():
    meridian_dist = load_meridian_distances()
    us32_dist = load_us32_distances()
    meridian_src, meridian_recv = load_meridian_coords()
    us32_src, us32_recv = load_us32_coords()
    print('Ex05 SKS record sections (radial vs transverse, jointly normalized):')
    print(f'  INPUT_ROOT = {INPUT_ROOT}')
    any_done = False
    for label, case_rel, geom, depth, tshift in SCENARIOS:
        stem = os.path.basename(case_rel).replace('sim_', '').replace('sim', '').strip('_')
        station_dir = resolve_station_dir(case_rel)
        n2d = us32_dist if geom == 'us32' else meridian_dist
        if plot_scenario(label, stem, station_dir, n2d, src_depth_km=depth,
                         right_shift_s=tshift, large_fonts=(geom == 'us32'),
                         filt_pmin=(12.0 if geom == 'meridian' else None)):
            print(f'  done: {label}')
            any_done = True
            if geom == 'us32':
                plot_geometry(label, stem, us32_src, us32_recv)
            else:
                plot_geometry(label, stem, meridian_src, meridian_recv)
    if not any_done:
        print('  no scenario data available yet — run the cases first (see README).')

    # SKS shear-wave-splitting "sticks" map (real Long et al. 2009 vs synthetic US32)
    try:
        import make_splitting_figure
        make_splitting_figure.main()
    except Exception as e:
        print(f'  [skip] splitting-sticks map: {e}')

    # Publication composite: transverse splitting signature + sticks map
    try:
        import make_composite_figure
        make_composite_figure.main()
    except Exception as e:
        print(f'  [skip] composite figure: {e}')

    print('Ex05 figures complete.')


if __name__ == '__main__':
    main()
