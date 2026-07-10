# Example 00: Global 1D (PREM)

**Author:** Jonathan Wolf (wolf@ucsc.edu)

## Overview

Global wave propagation for the 2011 Mw 5.8 Virginia earthquake (event
`VIRGINIA_201108231751A`) in the 1D anisotropic PREM reference model. The
synthetics are recorded at the Global Seismic Network (GSN) and at the USArray
Transportable Array (TA).

## Simulation

- **Dominant period:** 50 s (the coarse mesh shipped with this example).
- **Mesh:** `input/global_mesh__prem_ani__50s.e` (Exodus). It is pre-generated
  and committed, so no meshing step is needed to reproduce. To regenerate it
  with Salvus Mesher Lite:

  ```bash
  python -m salvus_mesh_lite.interface AxiSEM \
      --basic.model prem_ani --basic.period 50 \
      --output_file input/global_mesh__prem_ani__50s.e
  ```

- **Azimuthal resolution Nr:** `CONSTANT 5` (`type_Nr: CONSTANT`,
  `constant: 5`), the recommended value for a single axial moment-tensor source.
- **Source:** single moment tensor at `lat/lon = 37.91, -77.93`, depth 12 km
  (below the solid surface), moment tensor (CMTSOLUTION ZRT order)
  `[4.71e24, 3.81e22, -4.74e24, 3.99e23, -8.05e23, -1.23e24]` with `unit: 1e-7`
  (dyn·cm → N·m). The source-time function is a `GaussianSTF` with
  `half_duration: 50.` and `use_derivative_integral: ERF`.
- Record length 1800 s, `NEWMARK` integrator, `Courant_number: 1.0`.

All input parameters live in `input/` (`inparam.model.yaml`,
`inparam.nr.yaml`, `inparam.source.yaml`, `inparam.output.yaml`,
`inparam.advanced.yaml`). They reference the mesh and station files by relative
name only.

## Output

Defined in `input/inparam.output.yaml`: two station groups, both in the **RTZ**
(radial, transverse, vertical) coordinate frame, solid medium.

- **`global_seismic_network_GSN`** — GSN stations (`input/GSN.txt`), channels
  `[U]` (= U1, U2, U3 displacement), sampled at `DT`, format `ASCII_STATION`:
  one ascii file **per station** (`<NET>.<NAME>.ascii`) plus a shared
  `data_time.ascii`.
- **`USArray_transportable`** — USArray TA stations (`input/US_TA.txt`),
  channels `[U3, E_I1, R3]` (vertical displacement, strain trace, vertical
  curl/rotation), sampled every 5 s, format `ASCII_CHANNEL`: per-MPI-rank
  subfolders `dir_rank<N>/` each holding one file **per channel**
  (`U3.ascii`, `E_I1.ascii`, `R3.ascii`) and a `list_station.info`, plus a
  top-level `rank_station.info` mapping
  `MPI_RANK STATION_KEY STATION_INDEX_IN_RANK`.

Both groups are written under `output/stations/`.

## Reproduce

### 1. Mesh

Already provided as `input/global_mesh__prem_ani__50s.e`. See the **Simulation**
section above to regenerate it.

### 2. Large data

None. Everything needed to run is committed in `input/`.

### 3. Run

Build AxiSEM3D first and make the solver available as the `axisem3d` binary on
your `PATH`, or copy it into this folder as `./axisem3d`. Then, from this
example directory:

```bash
mpirun -np 4 axisem3d --input input --output output
```

Output is written to `output/` inside this folder. The shipped coarse 50 s mesh
runs comfortably on 4 MPI ranks; increase `-np` if you have more cores. On a
cluster, wrap the `mpirun` command in your scheduler's job script and load your
compiler, MPI, netCDF and HDF5 modules.

### 4. Figures

After the run, generate the publication figures from this example's `output/`:

```bash
python make_figures.py
```

`make_figures.py` reads the raw AxiSEM3D output under `output/` (set by the
`INPUT` variable at the top of the script) and the event/station coordinates
from `input/`. It imports the shared `pub_style.py` helper from the examples
root (colorblind-safe Okabe–Ito palette). Figures are written to `figures/` as
both PDF and PNG:

- **`ex00_station_map`** — global Robinson-projection map (Cartopy, coastlines
  only) showing all GSN stations as triangles and the Virginia 2011 epicenter as
  a star.
- **`ex00_seismograms_ANMO`** — 3-component (R/T/Z) displacement seismogram (in
  metres) at station `IU.ANMO`, titled with the great-circle epicentral
  distance Δ (in degrees).
- **`ex00_USArray_U3_timelapse`** — 2×4 time-lapse snapshot grid of the USArray
  **U3** (vertical displacement) wavefield over the array, snapshots evenly
  spaced across the record, `RdBu_r` color scale.
- **`ex00_USArray_R3_timelapse`** — same time-lapse grid for the **R3**
  (vertical rotation) wavefield.

If Cartopy is not installed the map is skipped and the remaining figures still
build. The figures need `numpy`, `pyyaml`, `matplotlib`, and (for the map)
`cartopy`.

## Interactive post-processing

`post_processing.ipynb` is an optional companion notebook that reads the same
`output/` and walks through the event map, GSN seismograms, ObsPy filtering, and
a USArray wavefield animation. Run the simulation first (step 3), then open the
notebook.
