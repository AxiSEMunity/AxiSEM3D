# Example 04: Simple 3-D shapes (blob / mantle plume)

This example shows how to superimpose simple 3-D heterogeneities (a spherical
blob, a plume-like cylinder + ellipsoid, and a field of randomly placed
spheres) on a 1-D background model in AxiSEM3D. Each shape is built in Python,
written to a `StructuredGridV3D` NetCDF model, and injected on top of a
homogeneous Exodus mesh. The directory holds three self-contained cases:
`example_input_cartesian` (a single Cartesian blob), `example_single_plume`
(a plume in a global domain), and `example_release_paper` (the wave-scattering
case from the AxiSEM3D release paper, with ~1000 random spheres). A
step-by-step tutorial, including how to inspect the 3-D models in ParaView,
is given in the notebook `Example_4_README.ipynb`.

## Contents
- `Example_4_README.ipynb` — step-by-step tutorial for all three cases.
- `src/` — Python library used to build the 3-D shapes (`model.py`,
  `injector.py`, `ellipsoid.py`, `cylinder.py`, `sphere.py`, `slab.py`,
  `object.py`, `global*.py`, `clrs.py`, `gen_scripts.py`).
- `images/` — tutorial screenshots (domain diagram, ParaView/NetCDF walkthrough).
- `figures/` — pre-rendered example figures.
- `make_figures.py` — regenerates the per-case figures directly from the
  simulation output under each case's `output/` directory.
- `make_blob_plume_figure.py` — publication composite (blob + plume wavefield
  time-lapses with their input velocity models); reads reduced cross-section
  archives (see [Figures](#figures)).
- `example_input_cartesian/` — single Cartesian blob:
  - `homogenous_cart.bm`, `gen_mesh.sh` — background model + mesh command.
  - `homogenous_cartesian.e` — shipped coarse Exodus mesh (period 0.2 s).
  - `gen_blob.py` — builds the 3-D sphere (`example_sphere_cartesian.nc`).
  - `inparam.*.yaml` — AxiSEM3D input (model / nr / source / output / advanced).
  - `stn.txt` — receiver list.
- `example_single_plume/` — plume in a global domain:
  - `homogenous_global.bm`, `gen_mesh.sh` — background model + mesh command.
  - `homogenous_global.e` — shipped coarse Exodus mesh (period 14 s).
  - `gen_plume.py` — builds the plume (`plume_1.nc`).
  - `gen_movie.py` — makes a ParaView movie from the wavefield output.
  - `inparam.*.yaml` — AxiSEM3D input.
- `example_release_paper/` — release-paper wave-scattering case:
  - `input/homogenous.bm`, `input/gen_mesh.sh` — background model + mesh command.
  - `input/paper_example.e` — shipped coarse Exodus mesh (period 0.08 s).
  - `randomblobs` — fixed list of sphere positions/sizes used in the paper.
  - `generate_3D_model.py` — builds the 3-D model (`paper_example.nc`) from
    `randomblobs`.
  - `input/inparam.*.yaml`, `input/paperstns.txt` — AxiSEM3D input + receivers.
  - `reproduce_paper_figure.py` — record-section figure from the output.

## Reproduce

These steps assume you have **built AxiSEM3D** and either have the `axisem3d`
binary on your `PATH` or have copied it into this example directory as
`./axisem3d` (see the top-level installation docs). The solver is always
launched with a portable MPI command of the form
`mpirun -np <N> axisem3d --input <case> --output <case>/output`.

### 1. Mesh
Each case ships its own coarse Exodus mesh (all well under 4 MB), so no mesh
needs to be regenerated to run the examples. The shipped meshes and the exact
commands that produced them (from each `gen_mesh.sh`) are:

- Cartesian blob — `example_input_cartesian/homogenous_cartesian.e`, period 0.2 s:
  ```bash
  python -m salvus_mesh_lite.interface AxiSEMCartesian \
    --basic.model homogenous_cart.bm --basic.period .2 \
    --cartesian2Daxisem.x 25. --cartesian2Daxisem.min_z 0.0 \
    --attenuation.frequencies 0.01 10. \
    --output_filename homogenous_cartesian.e
  ```
- Single plume — `example_single_plume/homogenous_global.e`, period 14 s:
  ```bash
  python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model homogenous_global.bm --basic.period 14 \
    --output_filename homogenous_global.e
  ```
- Release paper — `example_release_paper/input/paper_example.e`, period 0.08 s:
  ```bash
  python -m salvus_mesh_lite.interface AxiSEMCartesian \
    --basic.model homogenous.bm --basic.period .08 \
    --cartesian2Daxisem.x 30. --cartesian2Daxisem.min_z 0.0 \
    --attenuation.frequencies 0.01 16. \
    --output_filename paper_example.e
  ```

Mesh generation requires `salvus_mesh_lite` (see
https://axisem3d.readthedocs.io/en/latest/installation/mesher.html).

### 2. Large data
Each `inparam.model.yaml` references a 3-D NetCDF model that is **not** shipped
because it exceeds the 4 MB limit. Generate it locally with the Python script
in the corresponding case, or download it from the example gallery / Zenodo:

| Case | `nc_data_file` | Generate with | Approx. size |
|------|----------------|---------------|--------------|
| Cartesian blob | `example_sphere_cartesian.nc` | `python3 gen_blob.py` | ~70 MB (gallery) |
| Single plume | `plume_1.nc` | `python3 gen_plume.py` | ~71 MB (gallery) |
| Release paper | `paper_example.nc` | `python3 generate_3D_model.py` | ~4.9 GB (Zenodo) |

The release-paper model is large enough that it is hosted on Zenodo; the two
smaller models can be obtained from the gallery or regenerated in seconds with
the supplied scripts. `generate_3D_model.py` reads the shipped `randomblobs`
file so that the paper model is reproduced exactly.

### 3. Run
Once the 3-D NetCDF model exists in a case directory, point `--input` at that
case folder and launch the solver with `mpirun`. The shipped meshes are coarse,
so a handful of MPI ranks is plenty:

```bash
# (a) Cartesian blob  (~3-4 min on 4 cores)
cd example_input_cartesian
python3 gen_blob.py                              # -> example_sphere_cartesian.nc
mpirun -np 4 axisem3d --input . --output output
cd ..

# (b) Single global plume  (~13 min on 4 cores)
cd example_single_plume
python3 gen_plume.py                             # -> plume_1.nc
mpirun -np 4 axisem3d --input . --output output
cd ..
```

The **(c) release-paper** case is large (~5 GB model, ~60k steps) and is best
run on a cluster, but the launch command is identical — note that its inparam
files live in `input/`, so point `--input` at `example_release_paper/input`:

```bash
cd example_release_paper
python3 generate_3D_model.py                     # -> input/paper_example.nc (~5 GB)
mpirun -np 8 axisem3d --input input --output output
cd ..
```

On a cluster, wrap the `mpirun` command in your scheduler's job script and load
your compiler, MPI, netCDF and HDF5 modules.

### 4. Figures
Two plotting scripts live at the top of this example. Both import the shared
`pub_style` helper from `examples/pub_style.py` (located automatically, one
level up), run headless (matplotlib `Agg`), and write a PDF **and** a PNG into
`figures/`.

**`make_figures.py`** — reads each case's raw solver output directly from its
`output/` directory and regenerates the per-case figures. Run it after the
relevant case(s) have finished:

```bash
python3 make_figures.py
```

It produces, gated on the data being present:

| Figure (PDF + PNG) | Shows | Reads |
|--------------------|-------|-------|
| `figures/ex04_plume_cross_sections` | 2×2 vertical-displacement ($U_Z$) wavefield snapshots on the φ=0 slice for the plume | `example_single_plume/output/elements/...` |
| `figures/ex04_cartesian_seismograms` | RTZ station seismograms (m) for the Cartesian blob | `example_input_cartesian/output/stations/...` |
| `figures/ex04_release_paper_record_section` | Normalized vertical-displacement record section vs distance for the 1000 random spheres | `example_release_paper/output/stations/...` |

**`make_blob_plume_figure.py`** — the publication composite: a two-row figure
with the **(a)** Cartesian blob and **(b)** global plume wavefield time-lapses,
each preceded by its input velocity model, with the velocity anomalies, source
and free surface outlined.

```bash
python3 make_blob_plume_figure.py     # -> figures/ex04_blob_plume_snapshots.{pdf,png}
```

This script does **not** read the raw `output/` directly. It reads four reduced
NumPy archives from a `data/` directory next to the script (set by the
`INPUT_DIR` variable at the top of the file):

| Archive in `data/` | Content | Derived from |
|--------------------|---------|--------------|
| `cartesian_blob.npz` | blob input $V_S$ model (`x`, `z`, `dvs_xz`) | the `.nc` model from `gen_blob.py` |
| `plume_model.npz` | plume input $V$ model (`lat`, `depth`, `dv_xz`, …) | the `.nc` model from `gen_plume.py` |
| `cartesian_cross_section.npz` | blob $U_Z(s, z, t)$ on the φ=0 slice (`coords_sz`, `wave_zt`, `times`) | the blob run's `output/elements/` |
| `plume_cross_section.npz` | plume $U_Z(s, z, t)$ on the φ=0 slice | the plume run's `output/elements/` |

The two `*_cross_section.npz` files are **extracts** of the element NetCDF the
runs write under `output/elements/` (the raw element output is far larger than
the figure needs, so it is reduced to a compact archive first; the 30-km blob
slice is ~90 MB, well over the 4 MB repo limit, hence not shipped). To build
them, run cases (a) and (b) above, then extract the φ=0 ($U_Z$) slice from each
`output/elements/orthogonal_azimuthal_slices/` into the `coords_sz` / `wave_zt`
/ `times` arrays (the same reduction `make_figures.py` performs in-memory for
its plume panel). The two small model archives can be derived from the
generated `.nc` models. If any archive is missing, the script prints a `[skip]`
message and exits cleanly.

Other post-processing helpers:
- `example_single_plume/gen_movie.py` — builds a ParaView movie from the plume
  element output (reads `output/elements/orthogonal_azimuthal_slices`).
- `example_release_paper/reproduce_paper_figure.py` — the release-paper
  record-section figure; point its `ddir` at the run's station output
  (e.g. `output/stations/station_1`).
