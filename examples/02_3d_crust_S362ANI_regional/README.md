# Example 02: Regional 3-D crust (CRUST1.0) on S362ANI mantle

**Author:** Jonathan Wolf (wolf@ucsc.edu)

A regional simulation that compares a 1-D crust against a 3-D crust built from
CRUST1.0, both superimposed on the 3-D mantle model S362ANI over a 1-D PREM
(anisotropic) background. Running both cases shows how lateral crustal structure
reshapes regional waveforms across the USArray Transportable Array.

- Source: 2011 Virginia earthquake (Mw ~5.8), location [37.91, -77.93], depth 12 km
- Stations: USArray Transportable Array (`US_TA.txt`, 1834 stations)
- Mesh: regional chunk, CMB to surface, 60 degrees from the source, 50 s period

## Contents
- `input_share/` — inputs common to both cases: the coarse 50 s mesh
  `regional_mesh__prem_ani__50s.e`, `US_TA.txt` (receivers), and the
  `inparam.{advanced,nr,output,source}.yaml` configs.
- `input_with_1d_crust/inparam.model.yaml` — model config for the 1-D crust case
  (S362ANI mantle only).
- `input_with_3d_crust/inparam.model.yaml` — model config for the 3-D crust case;
  ships `crust1_geometric.nc` and `crust1_oceanload.nc` (the volumetric crust file
  is large, see **Large data**).
- `crust1.0_raw_data/` — raw CRUST1.0 grids (`crust1.bnds/.rho/.vp/.vs`) and their
  `readme.txt`, the source data for `gen_crust1.ipynb`.
- `gen_crust1.ipynb` — builds the CRUST1.0 NetCDF files in `input_with_3d_crust/`
  from `crust1.0_raw_data/`.
- `post_processing.ipynb`, `make_figures.py` — visualize and compare the results.
- `figures/` — pre-rendered figures (see **Reproduce → Figures**).

## Output
Defined in `input_share/inparam.output.yaml` (identical for both runs): one
station group `USArray_transportable` in the **RTZ** (radial, transverse,
vertical) frame, solid medium, channels `[U]` (= U1, U2, U3 displacement), the
1834 USArray TA stations from `US_TA.txt`, format `ASCII_CHANNEL`. That format
writes, under `output/stations/USArray_transportable/`, per-rank subfolders
`dir_rank<N>/` each holding one file per channel (`U1.ascii`, `U2.ascii`,
`U3.ascii`) and a `list_station.info`, plus a top-level `rank_station.info`
(mapping `MPI_RANK STATION_KEY ...`) and a shared `data_time.ascii`. The figures
use the **U3** (vertical displacement) channel.

## Reproduce

### 1. Mesh
The coarse mesh `input_share/regional_mesh__prem_ani__50s.e` (50 s period) is
shipped, so no meshing step is required. It spans CMB to surface vertically and
60 degrees from the source horizontally. To regenerate (or refine to a shorter
period by lowering `--basic.period`):
```
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani --basic.period 50 \
    --spherical.min_radius 3480 --chunk2D.max_colatitude 60 \
    --output_file regional_mesh__prem_ani__50s.e
```

### 2. Large data
Two model files referenced by the configs are too large to ship here and must be
obtained before running:
- `S362ANI_percent.nc` (~5 MB) — the 3-D mantle model, used by both cases.
  Download from IRIS EMC (model S362ANI) or the AxiSEM3D example gallery, and
  place it in `input_share/`.
- `crust1_volumetric.nc` (~30 MB) — the volumetric CRUST1.0 model used by the
  3-D crust case. Obtain it from the example gallery, or regenerate it with
  `gen_crust1.ipynb` from `crust1.0_raw_data/`, and place it in
  `input_with_3d_crust/`.

### 3. Run
Build AxiSEM3D and make the `axisem3d` binary available on your `PATH`, or copy
the built binary into this folder as `./axisem3d` (and call it as `./axisem3d`
below). Make sure your MPI / netCDF / HDF5 runtime is available.

Each case runs from a self-contained `input/` directory that combines the shared
inputs (`input_share/`) with that case's `inparam.model.yaml`. First place
`S362ANI_percent.nc` in `input_share/` as described in **Large data** (both
cases read it). Then assemble and run the two cases as follows, from this
example's directory.

**1-D crust case** (S362ANI mantle, 1-D PREM crust):
```
mkdir -p simu_with_1d_crust/input
cp input_share/*.yaml input_share/US_TA.txt \
   input_share/regional_mesh__prem_ani__50s.e input_share/S362ANI_percent.nc \
   simu_with_1d_crust/input/
cp input_with_1d_crust/inparam.model.yaml simu_with_1d_crust/input/

cd simu_with_1d_crust
mpirun -np 4 axisem3d --input input --output output
cd ..
```

**3-D crust case** (S362ANI mantle, 3-D CRUST1.0 crust):
```
mkdir -p simu_with_3d_crust/input
cp input_share/*.yaml input_share/US_TA.txt \
   input_share/regional_mesh__prem_ani__50s.e input_share/S362ANI_percent.nc \
   simu_with_3d_crust/input/
cp input_with_3d_crust/inparam.model.yaml simu_with_3d_crust/input/
# the 3-D crust NetCDF data the model config reads (crust1_volumetric.nc is the
# large file from **Large data**; place it in input_with_3d_crust/ first):
cp input_with_3d_crust/crust1_geometric.nc input_with_3d_crust/crust1_oceanload.nc \
   input_with_3d_crust/crust1_volumetric.nc simu_with_3d_crust/input/

cd simu_with_3d_crust
mpirun -np 4 axisem3d --input input --output output
cd ..
```

Results are written under `simu_with_<case>_crust/output/`. `-np 4` is fine for
the shipped coarse 50 s mesh on a laptop or workstation; raise it on a larger
machine or cluster. The 3-D crust run is the heavier of the two. On a cluster,
wrap the `mpirun` command in your scheduler's job script and load your compiler,
MPI, netCDF and HDF5 modules.

### 4. Figures
After both runs finish, render the comparison figures:
```
python make_figures.py
```
`make_figures.py` reads the two runs' station output directly from
`simu_with_1d_crust/output` and `simu_with_3d_crust/output` (the `INPUT_1D` /
`INPUT_3D` variables at the top of the script), parses the `ASCII_CHANNEL` U3
data, and writes PDF + PNG of each figure into `figures/`. It imports the shared
`pub_style.py` helper from the examples root (colorblind-safe Okabe–Ito palette,
metres for displacement) — no copy is needed, the script adds the examples root
to `sys.path` itself. If the simulation output is not present, the script prints
where it expected it and exits without error.

Figures written to `figures/`:
- `ex02_1D_vs_3D_crust_record_section.{pdf,png}` — vertical-displacement (U3)
  record section plotted against epicentral distance, overlaying the 1D-crust
  (blue) and 3D-crust (vermillion) traces at 60 stations spread evenly across
  the distance range. Shows how the 3-D crust reshapes and delays the regional
  wavetrain.
- `ex02_USArray_1D_crust_timelapse.{pdf,png}` — 2×4 time-lapse grid of the
  **1D-crust** USArray U3 wavefield over the array.
- `ex02_USArray_3D_crust_timelapse.{pdf,png}` — the same time-lapse grid for the
  **3D-crust** run. Both grids share one colour scale so the two are directly
  comparable side by side.
- `ex02_source_receiver_map.{pdf,png}` — orthographic map of the Virginia 2011
  source (star) and the USArray TA receivers (triangles) with great-circle
  paths. Requires Cartopy; if Cartopy is absent the script skips this panel and
  continues.

`post_processing.ipynb` is an alternative, interactive walk-through of the same
output (seismogram comparison plus optional wavefield animations).
