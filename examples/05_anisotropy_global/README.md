# Example 05: Full (3-D) anisotropy

This example demonstrates global seismic-wave propagation through fully
anisotropic Earth models, where the elastic response is described by the full
3-D stiffness tensor (`C11 ... C66`) rather than by transverse-isotropy
parameters alone. Several complementary workflows are collected here:

- **PREM anisotropy, two equivalent routes** — running an *anisotropic* PREM
  mesh directly versus running an *isotropic* PREM mesh and adding the same
  anisotropy back as a 3-D overlay. The two should reproduce one another.
- **Deep-mantle anisotropy** — an isotropic PREM mesh with a full-`Cij`
  anisotropic layer inserted in the lowermost mantle.
- **2012-07-03 paper example (US32)** — olivine (A-type) anisotropy beneath the
  USArray Transportable Array for the 2012 event, used to model SKS/SKKS
  shear-wave splitting. Two Nr (azimuthal-order) treatments are provided.

All shipped configurations use a coarse **50 s** mesh so the cases run quickly
on a workstation; the same inputs scale to higher-resolution meshes (see
**Mesh** below).

## Contents

```
README.md                         this file
make_figures.py                   SKS/SKKS record sections + geometry maps
make_splitting_figure.py          fast-axis splitting "sticks" map
make_composite_figure.py          publication composite (record section + map)
figures/                          reference figures produced by the scripts above
data/                             small metadata read by the plotting scripts
  sks_stations.txt                meridian receiver line (lat/lon)
  us32_stations.txt              US32 Transportable-Array receivers (lat/lon)
  splitting/                      per-station SKS-splitting measurements (phi, dt)

PREM_anisotropy_w_and_wo_full_Cij_50s/
  sim1_ani_prem_mesh/input/       anisotropic PREM mesh, run directly
  sim2_iso_prem_mesh_plus_ani/    isotropic PREM mesh + anisotropy overlay
                                  (prem_ani_upper.nc)
deep_mantle_anisotropy_full_Cij_50s/
  sim_lowermost_mantle_ani/input/ isotropic mesh + lowermost-mantle anisotropy
                                  (lowermost_mantle_ani.nc)
2012-07-03_paper_example_50s/
  sim_US32_olivineE/input/        US32 olivine anisotropy, POINTWISE Nr
  sim_US32_olivineE_fullNU/input/ same model, CONSTANT Nr variant
```

Each `input/` directory holds the standard AxiSEM3D parameter set
(`inparam.model.yaml`, `inparam.nr.yaml`, `inparam.source.yaml`,
`inparam.output.yaml`, `inparam.advanced.yaml`), the coarse Exodus mesh, the
receiver list, and any small 3-D-model NetCDF file the case needs. Receivers are
written as one station group (`JW_stations`) in `ASCII_STATION` format.

| Case | mesh | receiver file | 3-D model file | Nr |
| --- | --- | --- | --- | --- |
| sim1_ani_prem_mesh | `AxiSEM_prem_ani_50.e` | `stations.txt` | (in mesh) | CONSTANT |
| sim2_iso_prem_mesh_plus_ani | `AxiSEM_prem_iso_50.e` | `stations.txt` | `prem_ani_upper.nc` (shipped) | CONSTANT |
| sim_lowermost_mantle_ani | `AxiSEM_prem_iso_50.e` | `stations.txt` | `lowermost_mantle_ani.nc` (shipped) | CONSTANT |
| sim_US32_olivineE | `AxiSEM_prem_iso_50.e` | `eq3.txt` | `jons_azi_410_3.nc` (large) | POINTWISE (`scanning_output_Nr.nc`, large) |
| sim_US32_olivineE_fullNU | `AxiSEM_prem_iso_50.e` | `eq3.txt` | `jons_azi_410_3.nc` (large) | CONSTANT |

The `stations.txt` receiver grids were generated with the small
`make_stations.py` helper kept beside them.

## Mesh

Each case ships its own coarse **50 s** Exodus mesh in its `input/` directory:

- `AxiSEM_prem_iso_50.e` — isotropic PREM (6400 elements, ~0.4 MB)
- `AxiSEM_prem_ani_50.e` — radially anisotropic PREM (4176 elements, ~0.4 MB)

The meshes are built with the AxiSEM mesher (`salvus_mesh_lite`) from the
predefined PREM background models (`prem_iso`, `prem_ani`):

```bash
# isotropic PREM, 50 s -> 6400 elements
python3 -m salvus_mesh_lite.interface AxiSEM \
  --basic.model prem_iso --basic.period 50 \
  --output_filename AxiSEM_prem_iso_50.e --overwrite_file

# anisotropic PREM, 50 s -> 4176 elements
# (the prem_ani background model already carries the radial anisotropy;
#  elements_per_wavelength 1.5 reproduces the shipped element count)
python3 -m salvus_mesh_lite.interface AxiSEM \
  --basic.model prem_ani --basic.period 50 \
  --advanced.elements_per_wavelength 1.5 \
  --output_filename AxiSEM_prem_ani_50.e --overwrite_file
```

To refine, lower `--basic.period` (e.g. `--basic.period 5` for a 5 s mesh — the
resolution used for the publication figures). The resulting high-resolution
meshes are large (tens to hundreds of MB) and are not shipped here; obtain them
from the gallery or regenerate with the command above. The `exodus_mesh:` entry
in `inparam.model.yaml` must point at whichever mesh you use.

## Large data

These files exceed the in-repository size limit and are not shipped. Download
them from the example gallery (or regenerate them) and place them in the
relevant `input/` directory before running:

- `jons_azi_410_3.nc` (~19 MB) — full-`Cij` olivine-anisotropy model beneath
  the array, required by both US32 cases.
- `scanning_output_Nr.nc` (~18 MB) — pointwise Nr(s,z) field from a prior
  wavefield-scanning run, required only by `sim_US32_olivineE`
  (`type_Nr: POINTWISE`). Either fetch it from the gallery, or produce it by
  first running the case with `enable_scanning: true` and moving the resulting
  `output/scanning_output_Nr.nc` into `input/`. The `..._fullNU` variant uses
  `type_Nr: CONSTANT` and does not need this file.

The two PREM-overlay and deep-mantle cases need no external data: their 3-D
model files (`prem_ani_upper.nc`, `lowermost_mantle_ani.nc`) are small and ship
in the example.

## Reproduce

### 1. Mesh

Use the shipped coarse 50 s meshes as-is, or build a finer mesh as described in
**Mesh** above (the publication figures use a 5 s mesh).

### 2. Large data

Fetch `jons_azi_410_3.nc` (and, for `sim_US32_olivineE`, `scanning_output_Nr.nc`)
into the US32 `input/` directories as described in **Large data**. The PREM and
deep-mantle cases need nothing extra.

### 3. Run

Build AxiSEM3D first so the `axisem3d` solver binary exists. Put it on your
`PATH` (then call it as `axisem3d`), or copy the built binary into this example
directory and call it as `./axisem3d`. The commands below assume it is on your
`PATH`.

Run each case with a portable MPI launch. Each command reads that case's
`input/` directory and writes its results to a sibling `output/` directory next
to it; the plotting scripts later look for exactly these `output/` directories.
Four ranks are plenty for the coarse 50 s meshes:

```bash
# PREM anisotropy, run on the anisotropic mesh directly
mpirun -np 4 axisem3d \
  --input  PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh/input \
  --output PREM_anisotropy_w_and_wo_full_Cij_50s/sim1_ani_prem_mesh/output

# PREM anisotropy, isotropic mesh + anisotropy overlay (prem_ani_upper.nc)
mpirun -np 4 axisem3d \
  --input  PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani/input \
  --output PREM_anisotropy_w_and_wo_full_Cij_50s/sim2_iso_prem_mesh_plus_ani/output

# Deep-mantle anisotropy (lowermost_mantle_ani.nc)
mpirun -np 4 axisem3d \
  --input  deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani/input \
  --output deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani/output

# US32 olivine, pointwise Nr (needs jons_azi_410_3.nc + scanning_output_Nr.nc; see Large data)
mpirun -np 4 axisem3d \
  --input  2012-07-03_paper_example_50s/sim_US32_olivineE/input \
  --output 2012-07-03_paper_example_50s/sim_US32_olivineE/output

# US32 olivine, constant Nr (needs jons_azi_410_3.nc; see Large data)
mpirun -np 4 axisem3d \
  --input  2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU/input \
  --output 2012-07-03_paper_example_50s/sim_US32_olivineE_fullNU/output
```

On the 50 s meshes the PREM and deep-mantle cases run in minutes on a handful of
cores; the US32 cases are heavier because of the full-`Cij` olivine model. Raise
`-np` to match the cores you have available.

On a cluster, wrap the `mpirun` command in your scheduler's job script and load
your compiler, MPI, netCDF and HDF5 modules.

### 4. Figures

Run all three plotting scripts (each can also be run on its own):

```bash
python3 make_figures.py
```

The scripts import the shared `pub_style` helper from the parent `examples/`
directory and locate everything else **relative to themselves**. They read each
case's seismograms from its `output/stations/JW_stations/` directory (the
`ASCII_STATION` output you get by running the cases above); set the `EX05_INPUT`
environment variable if your runs live elsewhere. Each scenario with no readable
output is skipped automatically.

- **`make_figures.py`** — for every scenario, a 2-panel **Radial | Transverse**
  SKS/SKKS record section (`figures/ex05_SKS_<case>.{png,pdf}`) plus a
  source-receiver geometry map (`figures/ex05_geometry_<case>.{png,pdf}`).
  Radial and transverse are normalized jointly per station, so the transverse
  energy — the shear-wave-splitting signature of anisotropy — is shown to scale
  against the radial SKS pulse. This script also drives the two below.
- **`make_splitting_figure.py`** — the fast-axis splitting "sticks" map over the
  Pacific Northwest (`figures/ex05_splitting_sticks.{png,pdf}`): real Long et
  al. (2009) measurements (red) vs the AxiSEM3D US32 synthetic (blue), at the
  same length scale. Reads only the small `data/splitting/` files, so it works
  without running any case.
- **`make_composite_figure.py`** — the publication composite
  (`figures/ex05_composite_splitting.{png,pdf}`): panel (a) the US32 radial and
  transverse record section, panel (b) the splitting-sticks map. If a US32
  *isotropic* run is also present (drop it under
  `2012-07-03_paper_example_50s/sim_US32_isotropic/`), its flat transverse
  traces are drawn behind in black as the no-anisotropy reference.

Reference outputs are already in `figures/`.

**Note on resolution and output frame.** The reference figures in `figures/`
were produced from a fine **5 s** run whose station output was written in the
**RTZ** frame (`coordinate_frame: RTZ`, `channels: [U]`). The US32 cases shipped
here already use RTZ output and reproduce their figures directly. The meridian
(PREM / deep-mantle) cases ship with `coordinate_frame: ENZ`; to regenerate
their radial/transverse record sections set `coordinate_frame: RTZ` in their
`inparam.output.yaml` and refine the mesh to ~5 s before running. The plotting
script uses only the first three output columns (the `U` displacement vector) and
skips any case whose output it cannot interpret.

## Notes

Prepared/updated by Jonathan Wolf.
