# 01 Global 3D Simulation (S362ANI)

**Author:** Jonathan Wolf (wolf@ucsc.edu)

A global simulation using the 3D mantle model **S362ANI**, superimposed on a 1D
PREM (anisotropic) background. The 3D model file `S362ANI_percent.nc` was
downloaded from the IRIS EMC.

Basic information about this simulation:
- Source: 2011 Virginia earthquake (Mw 5.8)
- Stations: GSN global network + USArray transportable array
- Period: 50 s

## The solver binary

These instructions assume you have built AxiSEM3D and either have the `axisem3d`
binary on your `PATH`, or have copied it into this example directory as
`./axisem3d`. Replace `axisem3d` in the commands below with `./axisem3d` if you
copied it here. The launch command uses plain `mpirun`, so the example runs on
any Linux machine or cluster with an MPI installation; nothing here is tied to a
particular HPC system.

## Reproduce

### 1. Mesh

The coarse 50 s mesh `input/global_mesh__prem_ani__50s.e` is already provided, so
no meshing step is required. To regenerate it:

```bash
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model prem_ani --basic.period 50 \
    --output_file input/global_mesh__prem_ani__50s.e
```

The 3D structure is applied as a perturbation on top of this 1D mesh, not baked
into it (see `input/inparam.model.yaml`).

### 2. Large data

The 3D model file `S362ANI_percent.nc` (~4.9 MB) is referenced by
`input/inparam.model.yaml` but is too large to ship in this repository. Download
it from the **IRIS EMC** (model S362ANI), or obtain it from the example gallery,
and place it in `input/` before running.

### 3. Run

From this example directory:

```bash
mpirun -np 4 axisem3d --input input --output output
```

`-np 4` is a sensible default for the shipped coarse 50 s mesh; increase it on a
larger machine. The run takes roughly 35 minutes on 4 cores. Output is written to
`output/` inside this folder.

Station output is defined in `input/inparam.output.yaml` (two station groups,
both in the RTZ frame, solid medium):

- **`global_seismic_network_GSN`** — GSN stations (`input/GSN.txt`), channels
  `[U]`, format `ASCII_STATION`: one ascii file per station
  (`<NET>.<NAME>.ascii`) plus a shared `data_time.ascii`, under
  `output/stations/global_seismic_network_GSN/`.
- **`USArray_transportable`** — USArray TA stations (`input/US_TA.txt`), channels
  `[U3, E_I1, R3]`, sampled every 5 s, format `ASCII_CHANNEL`: per-MPI-rank
  subfolders `dir_rank<N>/` each with one file per channel, plus a top-level
  `rank_station.info`, under `output/stations/USArray_transportable/`.

*(Optional)* To enable the 1D-vs-3D comparison in the figures and notebook, also
run example `00_global_1D` so its `output/` exists.

### 4. Figures

```bash
python make_figures.py
```

`make_figures.py` reads this example's run output directly from `output/` (the
`INPUT` variable at the top of the script, set relative to the script) and writes
both PDF and PNG of each figure into `figures/`. It imports the shared
`pub_style.py` helper (colorblind-safe Okabe–Ito palette, seismograms in metres),
which is staged one directory up at `examples/pub_style.py`; the script adds that
location to `sys.path` relative to itself, so no installation or path editing is
needed. The script uses the headless Agg backend and runs without a display.

Figures produced:

- **`ex01_1D_vs_3D_ANMO`** — 3-component (R/T/Z) displacement seismogram at
  station `IU.ANMO`, overlaying the 3D S362ANI trace (solid, colored) against the
  1D PREM trace (dashed gray). The title shows the epicentral distance Δ
  (Virginia → Albuquerque). The 1D overlay is drawn only if
  `../00_global_1D/output/` is present; otherwise the 3D trace is plotted alone.
- **`ex01_USArray_U3_timelapse`** — time-lapse grid of the 3D USArray **U3**
  (vertical displacement) wavefield over the array, with a shared RdBu_r colour
  scale and one colorbar, showing the wave propagating across the transportable
  array.
- **`ex01_source_receiver_map`** — orthographic map of the source-receiver
  geometry: orange star = Virginia 2011 epicentre, blue triangles = USArray TA
  stations, with great-circle paths. Requires Cartopy; if Cartopy is not
  installed this figure is skipped (the others still render).

Reference copies of all three figures are in `figures/`.

## Interactive post-processing

`post_processing.ipynb` is an alternative, interactive walk-through that
visualizes the seismograms and builds a USArray wavefield animation, and (if
example `00_global_1D` has been run) overlays the 1D vs 3D comparison. It reads
the same `output/` directory.
