#!/usr/bin/env python3
"""
Example 04: Simple 3D Shapes — publication figures.
  Fig 1: Plume wavefield cross-section snapshots
  Fig 2: Cartesian blob station seismograms
"""
import sys, os
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from pub_style import apply_style, save_fig, CHANNEL_COLORS, CHANNEL_LABELS
import matplotlib.pyplot as plt

apply_style()

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, 'figures')


def read_element_output(data_dir):
    nc_fnames = sorted([f for f in os.listdir(data_dir) if 'axisem3d_synthetics.nc' in f])
    print(f'  Reading {len(nc_fnames)} file(s)')
    nc_files = [xr.open_dataset(os.path.join(data_dir, fn)) for fn in nc_fnames]
    na_grid = nc_files[0].data_vars['list_na_grid'].values.astype(int)
    data_time = nc_files[0].data_vars['data_time'].values
    xda_coords = xr.concat([f.data_vars['list_element_coords'] for f in nc_files], dim='dim_element')
    dict_data_wave = {}
    for nag in na_grid:
        dict_data_wave[nag] = xr.concat(
            [f.data_vars[f'data_wave__NaG={nag}'] for f in nc_files],
            dim=f'dim_element__NaG={nag}').values
    return na_grid, data_time, xda_coords.values, dict_data_wave


# ==== Figure 1: Plume cross-section snapshots ====
plume_elem_dir = os.path.join(SCRIPT_DIR, 'example_single_plume', 'output',
                              'elements', 'orthogonal_azimuthal_slices')
if os.path.isdir(plume_elem_dir):
    print('Ex04 Fig 1: Plume cross-section snapshots')
    na_grid, data_time, el_coords, data_wave = read_element_output(plume_elem_dir)
    nelem = el_coords.shape[0]
    ngll = el_coords.shape[1]
    ntime = len(data_time)
    phi_slices = np.radians(np.arange(0, 360, 90))
    nslice = len(phi_slices)
    coords_sz = el_coords.reshape((nelem * ngll, 2))
    wave_dim = 2  # Z
    islice = 0
    phi = phi_slices[islice]

    wave = np.ndarray((nelem * ngll, ntime))
    for ie in range(nelem):
        wave[ie * ngll:ie * ngll + ngll, :] = data_wave[nslice][ie, islice, :, wave_dim, :]

    s_km = coords_sz[:, 0] * np.cos(phi) / 1e3
    z_km = coords_sz[:, 1] / 1e3

    snapshots = [int(ntime * f) for f in [0.2, 0.4, 0.6, 0.8]]
    snapshots = [t for t in snapshots if t < ntime]
    vmax = np.percentile(np.abs(wave), 99.5)

    fig, axes = plt.subplots(2, 2, figsize=(7, 7))
    axes = axes.ravel()
    for i, tstep in enumerate(snapshots):
        ax = axes[i]
        sc = ax.scatter(s_km, z_km, s=0.3, c=wave[:, tstep],
                        vmin=-vmax, vmax=vmax, cmap='RdBu_r', rasterized=True)
        ax.set_title(f't = {data_time[tstep]:.1f} s', fontsize=10)
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Depth (km)')
        ax.set_aspect('equal')
    fig.suptitle('Global plume — Vertical displacement wavefield', fontsize=12, y=1.01)
    cbar = fig.colorbar(sc, ax=axes.tolist(), shrink=0.4, pad=0.02, aspect=30)
    cbar.set_label('$U_Z$ (m)', fontsize=9)
    save_fig(fig, os.path.join(FIG_DIR, 'ex04_plume_cross_sections'))
else:
    print('  Plume element output not found, skipping.')

# ==== Figure 2: Cartesian blob station seismograms ====
cart_dir = os.path.join(SCRIPT_DIR, 'example_input_cartesian', 'output', 'stations')
station_groups = ([d for d in os.listdir(cart_dir) if os.path.isdir(os.path.join(cart_dir, d))]
                  if os.path.isdir(cart_dir) else [])

if station_groups:
    print('Ex04 Fig 2: Cartesian blob seismograms')
    sg = station_groups[0]
    sg_dir = os.path.join(cart_dir, sg)
    time_file = os.path.join(sg_dir, 'data_time.ascii')
    time = np.loadtxt(time_file)
    st_files = sorted([f for f in os.listdir(sg_dir) if f.endswith('.ascii') and f != 'data_time.ascii'])

    nst = len(st_files)
    fig, axes = plt.subplots(3, 1, figsize=(7, 5), sharex=True)
    cmap = plt.cm.viridis
    for ist, sf in enumerate(st_files):
        d = np.loadtxt(os.path.join(sg_dir, sf))
        color = cmap(ist / max(1, nst - 1))
        label = sf.replace('.ascii', '')
        for ich, ch in enumerate('RTZ'):
            if ich < d.shape[1]:
                axes[ich].plot(time, d[:, ich], lw=0.6, color=color,
                               label=label if ich == 0 else None)

    for ich, ch in enumerate('RTZ'):
        axes[ich].set_ylabel(f'{CHANNEL_LABELS[ch]} (m)', fontsize=9)
        axes[ich].tick_params(direction='in')
        if ich < 2:
            axes[ich].tick_params(labelbottom=False)
    axes[0].set_title('Cartesian blob — Station seismograms', fontsize=11)
    axes[0].legend(loc='upper right', fontsize=6, ncol=2)
    axes[-1].set_xlabel('Time (s)')
    save_fig(fig, os.path.join(FIG_DIR, 'ex04_cartesian_seismograms'))
else:
    print('  Cartesian blob station output not found, skipping.')

# ==== Figure 3: Release paper - record section ====
rp_dir = os.path.join(SCRIPT_DIR, 'example_release_paper', 'output', 'stations')
if os.path.isdir(rp_dir):
    station_groups_rp = [d for d in os.listdir(rp_dir) if os.path.isdir(os.path.join(rp_dir, d))]
    if station_groups_rp:
        print('Ex04 Fig 3: Release paper record section')
        sg = station_groups_rp[0]
        sg_dir = os.path.join(rp_dir, sg)
        time = np.loadtxt(os.path.join(sg_dir, 'data_time.ascii'))

        stns = np.loadtxt(
            os.path.join(SCRIPT_DIR, 'example_release_paper', 'input', 'paperstns.txt'),
            dtype=str, skiprows=1)
        stn_names = [f'{r[1]}.{r[0]}' for r in stns]
        stn_dists = [float(r[2]) / 1e3 for r in stns]

        st_files = sorted([f for f in os.listdir(sg_dir) if f.endswith('.ascii')
                           and f != 'data_time.ascii' and 'rank_station' not in f])

        nst = len(st_files)
        if nst > 0:
            dists = []
            traces = []
            for sf in st_files:
                name = sf.replace('.ascii', '')
                d = np.loadtxt(os.path.join(sg_dir, sf))
                idx_match = [i for i, n in enumerate(stn_names) if n == name]
                dist = stn_dists[idx_match[0]] if idx_match else 0
                vert = d[:, 2] if d.ndim > 1 else d
                dists.append(dist)
                traces.append(vert)

            sort_idx = np.argsort(dists)
            fig, ax = plt.subplots(figsize=(7, 8))
            for i, idx in enumerate(sort_idx):
                norm_val = max(np.max(np.abs(traces[idx])), 1e-30)
                ax.plot(time, traces[idx] / norm_val + i * 2,
                        color=CHANNEL_COLORS['Z'], lw=0.4)
                if i % 5 == 0:
                    ax.text(time[0] - 20, i * 2, f'{dists[idx]:.0f} km',
                            ha='right', va='center', fontsize=5)
            ax.set_xlabel('Time (s)')
            ax.set_yticks([])
            ax.set_ylim(-2, len(sort_idx) * 2)
            ax.spines['left'].set_visible(False)
            ax.set_title('1000 random spheres — Vertical displacement', fontsize=11)
            save_fig(fig, os.path.join(FIG_DIR, 'ex04_release_paper_record_section'))

print('Ex04 figures complete.')
