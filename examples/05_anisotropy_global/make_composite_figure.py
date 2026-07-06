#!/usr/bin/env python3
"""
Example 05: Global anisotropy — publication composite (US32, 2012-07-03 NZ event).

A single well-ordered figure with all legends OUTSIDE the axes:

  (a) Transverse-component record section. The ISOTROPIC (no-anisotropy) transverse
      traces are drawn behind, in black; the ANISOTROPIC olivine transverse traces are
      overlaid in blue. A purely isotropic Earth has zero transverse SKS, so the black
      traces are ~flat and the blue wiggles are the shear-wave-splitting signature. Both
      are normalized per station by the olivine radial SKS peak (so the split energy is
      shown to scale), and the synthetics are shifted +5 s so the SKS onset sits on the
      predicted line. The isotropic overlay is drawn only if a US32 isotropic run is
      available (see ISO_REL below); otherwise just the olivine transverse is shown.

  (b) The fast-axis splitting "sticks" map (real Long et al. 2009 vs AxiSEM3D synthetic),
      drawn by make_splitting_figure.draw_map.

INPUT: the US32 olivine (and, optionally, US32 isotropic) ASCII_STATION seismograms,
resolved relative to this example (run output/, or a flat data/ extraction), plus the
small splitting metadata in data/splitting/. Standalone: run
``python3 make_composite_figure.py`` from this folder.
"""
import os
import sys
import glob

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # examples/ root holds pub_style.py

from pub_style import apply_style, save_fig, COLORS
import matplotlib.pyplot as plt
from obspy.taup import TauPyModel
import make_splitting_figure as msf

apply_style()
plt.rcParams['figure.constrained_layout.use'] = False

INPUT_ROOT = os.environ.get('EX05_INPUT', HERE)
DATA_DIR = os.path.join(HERE, 'data')
MODEL = TauPyModel('prem')

EVENT = (-40.023, 173.756)     # 2012-07-03 NZ event lat/lon
DEPTH = 229.8                  # km
WIN_PRE = 15.0
SHIFT = 5.0                    # +s: put SKS onset (not peak) on the predicted line

# case-relative paths to the anisotropic (olivine) and isotropic US32 runs. The example
# ships the olivine case; an isotropic US32 run is optional (drop it under the path below
# to enable the black "isotropic transverse only" reference).
OLIV_REL = '2012-07-03_paper_example_50s/sim_US32_olivineE'
ISO_REL = '2012-07-03_paper_example_50s/sim_US32_isotropic'

C_OLIV = COLORS['primary']        # blue  = AxiSEM3D anisotropic synthetic
C_ISO = 'black'                   # black = isotropic reference (behind)
FS = dict(title=16, axlabel=15, tick=13, legend=12.5, panel=18)


def gc_deg(la1, lo1, la2, lo2):
    a = map(np.radians, (la1, lo1, la2, lo2)); la1, lo1, la2, lo2 = a
    return np.degrees(np.arccos(np.clip(
        np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1), -1, 1)))


def resolve_station_dir(case_rel):
    """Find the directory with the <station>.ascii files for a case (run output/ first,
    then a flat data/ extraction). Returns None if nothing is found."""
    base = os.path.basename(case_rel)
    for c in (os.path.join(INPUT_ROOT, case_rel, 'output', 'stations', 'JW_stations'),
              os.path.join(INPUT_ROOT, 'data', base),
              os.path.join(INPUT_ROOT, base)):
        if os.path.isfile(os.path.join(c, 'data_time.ascii')):
            return c
    # any group name under output/stations/
    sd = os.path.join(INPUT_ROOT, case_rel, 'output', 'stations')
    if os.path.isdir(sd):
        for grp in sorted(glob.glob(os.path.join(sd, '*'))):
            if os.path.isfile(os.path.join(grp, 'data_time.ascii')):
                return grp
    return None


def load_distances():
    f = os.path.join(DATA_DIR, 'us32_stations.txt')
    n2d = {}
    for line in open(f):
        c = line.split()
        if len(c) >= 4 and not line.startswith('#'):
            d = gc_deg(*EVENT, float(c[2]), float(c[3]))
            n2d[f'{c[1]}.{c[0]}'] = d
            n2d[c[0]] = d
    return n2d


def phase_time(d, phase):
    arr = MODEL.get_travel_times(source_depth_in_km=DEPTH, distance_in_degree=d,
                                 phase_list=[phase])
    return arr[0].time if arr else None


def load_scenario(case_rel, n2d):
    station_dir = resolve_station_dir(case_rel)
    if station_dir is None:
        return None, {}
    time = np.loadtxt(os.path.join(station_dir, 'data_time.ascii'))
    recs = {}
    for f in sorted(glob.glob(os.path.join(station_dir, '*.ascii'))):
        nm = os.path.basename(f)[:-6]
        if nm == 'data_time' or nm not in n2d:
            continue
        data = np.loadtxt(f)
        if data.ndim == 1 or data.shape[1] < 3:
            continue
        recs[nm] = (n2d[nm], data[:, :3])     # (dist, data[R,T,Z])
    return time, recs


def measure_offset(time, recs, half=12.0):
    offs = []
    for nm, (dist, data) in recs.items():
        t_sks = phase_time(dist, 'SKS')
        if t_sks is None:
            continue
        w = (time >= t_sks - half) & (time <= t_sks + half)
        if not w.any():
            continue
        env = np.maximum(np.abs(data[:, 0][w]), np.abs(data[:, 1][w]))
        if env.max() <= 0:
            continue
        offs.append((time[w] - t_sks)[int(np.argmax(env))])
    return float(np.median(offs)) if offs else 0.0


def panel_records(axR, axT):
    """Left axR: RADIAL component, anisotropic (olivine). Right axT: TRANSVERSE, isotropic
    (black, behind, if available) vs anisotropic olivine (blue). Both components are
    normalized per station by the SAME olivine radial SKS peak. Returns legend handles."""
    n2d = load_distances()
    t_o, rec_o = load_scenario(OLIV_REL, n2d)
    t_i, rec_i = load_scenario(ISO_REL, n2d)
    have_iso = bool(rec_i)
    common = sorted(rec_o, key=lambda nm: rec_o[nm][0])
    if have_iso:
        common = sorted(set(rec_o) & set(rec_i), key=lambda nm: rec_o[nm][0])
    if len(common) > 32:
        idx = np.linspace(0, len(common) - 1, 32).astype(int)
        common = [common[i] for i in idx]
    dists = [rec_o[nm][0] for nm in common]
    span = max(max(dists) - min(dists), 1.0)
    gain = span / 22.0

    offset = measure_offset(t_o, rec_o)
    far = max(dists)
    ts_far, tk_far = phase_time(far, 'SKS'), phase_time(far, 'SKKS')
    win_post = 45.0
    if ts_far and tk_far:
        win_post = float(min(max(45.0, (tk_far - ts_far) + 12.0), 165.0))

    skks_pts = []
    for nm in common:
        dist = rec_o[nm][0]
        t_sks = phase_time(dist, 'SKS')
        if t_sks is None:
            continue
        t_align = t_sks + offset
        w = (t_o >= t_align - WIN_PRE - SHIFT) & (t_o <= t_align + win_post)
        if not w.any():
            continue
        tw = t_o[w] - t_align + SHIFT
        Ro = rec_o[nm][1][:, 0][w]
        To = rec_o[nm][1][:, 1][w]
        norm = max(np.max(np.abs(Ro)), 1e-30)            # olivine radial SKS peak
        # RADIAL — anisotropic olivine only
        axR.plot(tw, Ro / norm * gain + dist, color=C_OLIV, lw=0.9, zorder=3)
        # TRANSVERSE — isotropic (behind, if available) then olivine, same time base & norm
        if have_iso:
            wi = (t_i >= t_align - WIN_PRE - SHIFT) & (t_i <= t_align + win_post)
            Ti = rec_i[nm][1][:, 1][wi]
            twi = t_i[wi] - t_align + SHIFT
            axT.plot(twi, Ti / norm * gain + dist, color=C_ISO, lw=1.3, zorder=2)
        axT.plot(tw, To / norm * gain + dist, color=C_OLIV, lw=0.9, zorder=3)
        t_skks = phase_time(dist, 'SKKS')
        if t_skks is not None and -WIN_PRE <= (t_skks - t_sks) <= win_post:
            skks_pts.append((t_skks - t_sks, dist))

    sp = np.array(sorted(skks_pts, key=lambda p: p[1])) if skks_pts else None
    for ax, ttl in ((axR, 'Radial (anisotropic)'), (axT, 'Transverse')):
        ax.axvline(0.0, color=COLORS['gray'], lw=1.0, ls='--', zorder=1)
        if sp is not None:
            ax.plot(sp[:, 0], sp[:, 1], color=COLORS['quaternary'], lw=1.3, ls='--', zorder=1)
        ax.set_xlim(-WIN_PRE, win_post)
        ax.set_ylim(min(dists) - gain * 1.5, max(dists) + gain * 1.5)
        ax.set_title(ttl, fontsize=FS['title'])
        ax.tick_params(labelsize=FS['tick'])
    axR.set_ylabel('Epicentral distance ($^\\circ$)', fontsize=FS['axlabel'])
    axT.tick_params(labelleft=False)

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=C_OLIV, lw=2.6, label='Anisotropic (olivine)')]
    if have_iso:
        handles.append(Line2D([0], [0], color=C_ISO, lw=2.6,
                              label='Isotropic (transverse only)'))
    handles += [
        Line2D([0], [0], color=COLORS['gray'], lw=1.2, ls='--', label='SKS (predicted)'),
        Line2D([0], [0], color=COLORS['quaternary'], lw=1.4, ls='--', label='SKKS (predicted)'),
        Line2D([0], [0], color='none', label='synthetics +5 s (onset on line)'),
    ]
    return handles


def main():
    try:
        import cartopy.crs as ccrs
    except ImportError:
        print('  (cartopy missing, skipping composite figure)')
        return
    if resolve_station_dir(OLIV_REL) is None:
        print('  [skip] composite: US32 olivine data not found (run the case first)')
        return

    fig = plt.figure(figsize=(16.0, 8.8))

    # (a) two shared-y record sub-panels: Radial (olivine) | Transverse (oliv [+ iso])
    rect_R = [0.050, 0.275, 0.150, 0.600]
    rect_T = [0.205, 0.275, 0.150, 0.600]
    ax_R = fig.add_axes(rect_R)
    ax_T = fig.add_axes(rect_T, sharey=ax_R)
    handles_a = panel_records(ax_R, ax_T)

    # (b) splitting sticks map (reuses the standalone map drawing)
    proj = ccrs.Mercator(central_longitude=0.5 * (msf.LON_MIN + msf.LON_MAX),
                         min_latitude=msf.LAT_MIN, max_latitude=msf.LAT_MAX)
    rect_b = [0.470, 0.125, 0.400, 0.750]
    ax_b = fig.add_axes(rect_b, projection=proj)
    msf.draw_map(fig, ax_b, rect_b, reference=False)   # scale now lives in the legend below
    ax_b.set_title('Fast-axis splitting map', fontsize=FS['title'])

    fig.text(rect_b[0] + 0.5 * rect_b[2], 0.095, 'Longitude',
             ha='center', va='center', fontsize=FS['axlabel'])
    fig.text(rect_b[0] - 0.052, rect_b[1] + 0.5 * rect_b[3], 'Latitude',
             ha='center', va='center', rotation=90, fontsize=FS['axlabel'])

    a_cx = 0.5 * (rect_R[0] + rect_T[0] + rect_T[2])      # centre of the pair
    fig.text(a_cx, 0.215, 'Time relative to predicted SKS (s)',
             ha='center', va='center', fontsize=FS['axlabel'])

    # --- legends, both OUTSIDE the axes, in the bottom margin ---
    fig.legend(handles=handles_a, loc='upper center', bbox_to_anchor=(a_cx, 0.175),
               ncol=2, frameon=False, fontsize=FS['legend'], handlelength=1.8,
               columnspacing=1.6)
    yb = (0.072 - rect_b[1]) / rect_b[3]
    hl = msf.one_second_handlelength(16.0, rect_b[2], FS['legend'])
    leg_b = ax_b.legend(handles=msf.legend_handles(), loc='upper center',
                        bbox_to_anchor=(0.5, yb), ncol=2, frameon=False,
                        fontsize=FS['legend'], handlelength=hl, columnspacing=2.0,
                        handletextpad=0.6, labelspacing=0.3, borderpad=0.2,
                        title='Stick length = 1 s splitting')
    leg_b.get_title().set_fontsize(FS['legend'])

    fig.text(rect_R[0] - 0.038, rect_R[1] + rect_R[3] + 0.030, '(a)',
             fontsize=FS['panel'], fontweight='bold', va='bottom')
    fig.text(rect_b[0] - 0.028, rect_b[1] + rect_b[3] + 0.030, '(b)',
             fontsize=FS['panel'], fontweight='bold', va='bottom')

    fig.suptitle('SKS shear-wave splitting — 2012-07-03 New Zealand event across '
                 'the US32 array', fontsize=17, y=0.965)
    save_fig(fig, os.path.join(HERE, 'figures', 'ex05_composite_splitting'), tight=False)
    print('  done: ex05_composite_splitting')


if __name__ == '__main__':
    main()
