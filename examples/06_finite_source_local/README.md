# Example 06: Kinematic finite-fault rupture (San Francisco Bay Area)

This example simulates ground motion from an earthquake in the San Francisco Bay
Area (SFBA) using a local 3D crustal velocity model with topography and a Moho
interface on a Cartesian (`AxiSEMCartesian`) mesh. A kinematic finite-fault
rupture can be approximated as a collection of point moment-tensor sources (one
Gaussian source-time function per subsource); the shipped `inparam.source.yaml`
uses a single representative moment-tensor source, and `input_setup.ipynb` shows
how to format a multi-subsource list. Surface ground motion is recorded on a
dense grid of stations for animating the wavefield and mapping peak ground
velocity.

## Contents
- `input/inparam.model.yaml` — 1D/3D model: references the Exodus mesh and the
  SFBA 3D velocity layers, topography and Moho undulation.
- `input/inparam.source.yaml` — time axis (`record_length: 30 s`) and the
  moment-tensor source; a Gaussian STF with `half_duration: 2 s` band-limited to
  the mesh.
- `input/inparam.nr.yaml` — azimuthal resolution Nr (POINTWISE, from
  `point_source_approx.nc`), with wavefield scanning enabled.
- `input/inparam.output.yaml` — station group `sfba_local` (reads
  `STATIONS_OUTPUT.txt`).
- `input/inparam.advanced.yaml` — advanced/solver settings.
- `input/AxiSEMCartesian_sfba_2018_2s.e` — the shipped 2 s Exodus mesh (33,900
  elements), matched to the source's `half_duration: 2 s` (a finer mesh only
  adds unresolved cost).
- `input/sfba_2018.bm` — polynomial background model used to build the mesh.
- `input/genmesh.sh` — command that generates the mesh from the `.bm` file.
- `input/STATIONS_OUTPUT.txt` — surface station grid for the recorded output.
- `input/example_input_source_list.txt` — example formatted subsource list
  produced by `input_setup.ipynb` (format template for a multi-subsource rupture).
- `input/point_source_approx.nc` — pointwise Nr field used by `inparam.nr.yaml`.
- `input/moho_smooth_24km.nc` — Moho undulation grid used by `inparam.model.yaml`.
- `input_setup.ipynb` — builds the source list and `STATIONS_OUTPUT.txt`.
- `make_figures.py` — generates the publication figures from `output/`.
- `make_movie.py` — generates `movies/ex06_sfba_velocity_Z.mp4`, an animation of
  the vertical-velocity wavefield over the full record.
- `post_processing.ipynb` — an older notebook alternative for the figures/
  animation; prefer `make_figures.py`/`make_movie.py` (it needs `dask`, which
  the `axisem3d` conda env doesn't have, and used a stale nm→m unit conversion).

## Requirements

You need a built `axisem3d` solver and an MPI runtime. Either put the `axisem3d`
binary on your `PATH` or copy it into this directory as `./axisem3d`. The figure
script additionally needs Python with `numpy`, `pandas`, `netcdf4`, `scipy`,
`matplotlib` and `cartopy` (no `xarray`/`dask` — the per-rank files are read
directly with `netCDF4` and concatenated by hand, since `dask` is often
unavailable on HPC Python environments).

## Reproduce

### 1. Mesh
The shipped mesh `input/AxiSEMCartesian_sfba_2018_2s.e` is a **2 s** Cartesian
`AxiSEMCartesian` mesh (33,900 elements), built from the polynomial background
model `input/sfba_2018.bm`. It is already included, so you only need this step
to regenerate or refine it (lower `--basic.period` for a finer mesh — but note
the source's `half_duration` is 2 s, so a finer mesh resolves nothing extra and
only adds cost; a 0.5 s mesh of this same model runs to ~475k elements and is
also more prone to boundary-reflection instability at the default Courant
number). Run from `input/`:
```
sh genmesh.sh
```
which executes:
```
python -m salvus_mesh_lite.interface AxiSEMCartesian \
    --basic.model sfba_2018.bm --basic.period 2.0 \
    --attenuation.frequencies 0.01 10.0 \
    --cartesian2Daxisem.x 50.0 \
    --output_file AxiSEMCartesian_sfba_2018_2s.e
```

### 2. Large data
The 3D velocity model referenced in `inparam.model.yaml` is too large to ship in
git and must be downloaded and copied into `input/`:
`SFBA_layer1_min500.nc`, `SFBA_layer2.nc`, `SFBA_layer3.nc`, `SFBA_layer4.nc`,
and `topo_smooth_100m_upsampled.nc`. Obtain them from Zenodo:
https://zenodo.org/records/20017355

### 3. Run
The shipped 2 s mesh has 33,900 elements. From this example directory, run the
solver inline:
```
mpirun -np 32 axisem3d --input input --output output
```
(32 ranks leaves ~1,000 elements/rank — fine for a workstation; scale the rank
count up on a cluster, e.g. one rank per ~100–150 elements.) On a cluster, wrap
the `mpirun`/`srun` command in your scheduler's job script and load your
compiler, MPI, netCDF and HDF5 modules. The shipped `Courant_number: 0.5` in
`inparam.source.yaml` gives a safe stability margin for this mesh/model; a
value of 0.6 was found to be marginally unstable (a late-onset, growing
instability once the wavefront reaches the ±50 km domain boundary, around
t ≈ 8 s) — don't raise it without retesting.

The run writes per-rank NetCDF station output to
`output/stations/sfba_local/axisem3d_synthetics.nc.rank*`.

### 4. Figures
Once `output/` exists, generate the publication figures:
```
python make_figures.py
```
`make_figures.py` reads the per-rank station NetCDF directly from
`output/stations/sfba_local/` (the `INPUT` variable at the top of the script,
defaulting to `./output`) and `input/STATIONS_OUTPUT.txt`. It differentiates the
recorded displacement in time to velocity, then writes three figures (PDF + PNG)
to `figures/`:

- `figures/ex06_wavefield_snapshots.{pdf,png}` — 4-panel map snapshot at
  t ≈ 60 % of the record: vertical (Z), east (E) and north (N) velocity
  (RdBu_r, symmetric scale at the data's 99.5th percentile) plus a running
  peak-ground-velocity panel (YlOrRd, its own scale from the PGV field's 99th
  percentile — PGV is a running max, so it saturates a shared velocity scale),
  in a UTM zone 11 projection over the SFBA. Each field is linearly
  interpolated from the station grid onto a continuous surface that is only
  drawn inside the stations' convex hull (i.e. exactly where the simulation
  actually sampled it — no extrapolation beyond the recorded disk).
- `figures/ex06_pgv_map.{pdf,png}` — peak ground velocity (max of the velocity
  norm over the full record) at every station, same continuous-field
  treatment and PGV-specific scale, YlOrRd, same SFBA extent.
- `figures/ex06_velocity_seismograms.{pdf,png}` — three-component (E/N/Z)
  velocity seismograms at 5 high-PGV stations, each labeled with its (x, y)
  location in km.
- `movies/ex06_sfba_velocity_Z.mp4` — vertical-velocity wavefield animation
  (`make_movie.py`), same continuous-field/coastline treatment, all 176 frames
  of the record at 12 fps.

All figures embed **editable Helvetica** text (`pdf.use14corefonts`, core Type1
font — not an embedded subset — so text stays fully editable in Illustrator).

The maps use coastlines only (Cartopy `resolution='10m'`, no state or country
borders) on a white background — this local ~90 km-wide extent is zoomed in
enough that 10m detail (SF Bay's inlets/islands) is legible and worthwhile,
unlike a global or continental map, where a coarser 50m/110m resolution is the
better (faster, still-accurate-at-that-scale) choice. The figure script
imports the shared
`pub_style.py` helper from the parent `examples/` directory (colorblind-safe
palette, SI units); if the simulation output is absent it prints a notice and
exits cleanly. `post_processing.ipynb` is a notebook alternative that also builds
the surface-wavefield animation.

## Notes
The moment-tensor `data` values in `inparam.source.yaml` are the GCMT tensor in
**dyne·cm** — hence `unit: 1.0e-7` to convert to N·m (scalar M0 → M_w 4.4). Using
`unit: 1` instead reads the same numbers as already being in N·m, giving M_w 9.0:
a source ~1e7× too strong. If you swap in your own moment tensor, double check
which convention it's quoted in before setting `unit`.

Prepared/updated by Jonathan Wolf.
