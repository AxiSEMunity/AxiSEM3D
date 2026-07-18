# Example 11: Wave visualization on the Moon (meteoroid impact)

Forward AxiSEM3D simulation of the seismic surface wavefield excited by a shallow
meteoroid impact on the Moon, using the 1-D anelastic lunar model **ISSI_MOON_M1**
(Garcia et al. 2019). The source is an on-axis isotropic explosion (the volumetric
excavation component of an impact). Because the source sits on the symmetry axis,
the radiated wavefield is **axisymmetric** (azimuthal order m = 0): the entire
lunar surface field is captured by a single dense meridian line of receivers and
reconstructed everywhere in post-processing by rotating that line about the impact
point. This example demonstrates how the long, slowly-decaying coda characteristic
of the Moon's strongly scattering, low-attenuation interior is produced, and how
an axisymmetric run reproduces a global surface wavefield at a fraction of the cost
of a full 3-D simulation.

## Contents

```
11_wave_visualization_Moon/
├── README.md                       # this file
├── README_gallery_pyvista_viz.md   # upstream gallery doc for the optional 3-D PyVista visualization
├── make_figures.py                 # publication figures (record section, snapshots, seismogram)
├── conda_environment.yaml          # conda env for the optional PyVista visualization
├── combine_png.ipynb               # notebook helper for the PyVista movie frames
├── stations_processing.py          # PyVista pipeline part 1 (build wavefield meshes)
├── png_creation_seismo.py          # PyVista pipeline part 2 (render frames)
├── postprocessing_util.py          # helpers for the PyVista pipeline
├── postprocessing_util_observations.py
├── stations_processing_1D_3D.py
├── multiprocessing/                # parallel variants of the PyVista station/png scripts
├── figures/                        # reference output figures (PDF + PNG)
└── sim/
    ├── ISSI_MOON_M1_atten.bm       # 1-D anelastic lunar model (Garcia et al. 2019)
    └── input/                      # the shipped coarse run
        ├── moon_issi_m1_15s.e      # Exodus mesh, 15 s minimum period, 18,400 elements (0.8 MB)
        ├── moon_line.txt           # 1799 receivers, 0.1° spacing along the equator (lon = epicentral distance Δ)
        ├── inparam.model.yaml      # mesh file, geodesy (SPHERE), absorbing boundary, attenuation CG4
        ├── inparam.source.yaml     # isotropic explosion at 1 km depth, GaussianSTF
        ├── inparam.nr.yaml         # Nr = CONSTANT 5
        ├── inparam.output.yaml     # surface_line station group, RTZ U, NETCDF
        └── inparam.advanced.yaml   # verbosity / mpi / misc
```

The model file is `ISSI_MOON_M1_atten.bm`: a 1-D anelastic, isotropic lunar model
(radius R = 1737.1 km, crust + mantle and a fluid core; columns are
`depth vp vs rho qmu qkappa`, units km). The source is a `MOMENT_TENSOR` isotropic
explosion `data: [1e15, 1e15, 1e15, 0, 0, 0]` at 1 km depth below the solid surface
at `latitude_longitude: [0, 0]`, with a `GaussianSTF` source-time function and
`Courant_number: 0.5` (a single vertical force triggered a near-axis instability in
this lunar mesh; the isotropic explosion is stable and axisymmetric).

## Mesh

The shipped mesh `sim/input/moon_issi_m1_15s.e` is a **coarse 15 s** AxiSEM mesh
(18,400 elements, 0.8 MB) built from the `.bm` model with `salvus_mesh_lite`. To
regenerate it, run from `sim/` with a Python that has `salvus_mesh_lite`:

```bash
cd sim
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model ISSI_MOON_M1_atten.bm --basic.period 15 \
    --output_filename input/moon_issi_m1_15s.e
```

To **refine** the mesh, lower `--basic.period` (and lengthen the source-time
function `half_duration` in `inparam.source.yaml` to resolve the new minimum
period). The high-resolution 4 s / 2 s / 1 s meshes used for the headline figures
exceed the repository size cap and are listed under **Large data** below.

## Large data

The high-resolution Moon meshes are not shipped here (they exceed the 4 MB
per-file limit). Regenerate any of them with the `salvus_mesh_lite` command above
by changing `--basic.period`, or obtain the pre-built `.e` files from the gallery:

| Mesh | Period | Elements | Approx. size |
|------|--------|----------|--------------|
| `moon_issi_m1_4s.e` | 4 s  | 251,304   | ~9.8 MB  |
| `moon_issi_m1_2s.e` | 2 s  | ~1.0 M    | ~38 MB   |
| `moon_issi_m1_1s.e` | 1 s  | ~3.99 M   | ~158 MB  |

The raw per-rank seismogram output of a high-resolution run
(`output_*s/stations/surface_line/axisem3d_synthetics.nc.*`, several GB total)
belongs on Zenodo, not in the repository.

## Reproduce

### 1. Mesh

The shipped run uses the coarse 15 s mesh `sim/input/moon_issi_m1_15s.e`
(18,400 elements, 0.8 MB), which is already in the repository. Regenerate or
refine it with the `salvus_mesh_lite` command in the [Mesh](#mesh) section above.

### 2. Large data

Only needed for the high-resolution headline figures. The refined
4 s / 2 s / 1 s meshes and the multi-GB raw station output of those runs are
**not** shipped (see [Large data](#large-data)). Regenerate the meshes with
`salvus_mesh_lite`, or obtain the pre-built `.e` files and reduced outputs from
the AxiSEM3D gallery / Zenodo. The shipped coarse run needs none of this.

### 3. Run

Build AxiSEM3D and either put the `axisem3d` binary on your `PATH` or copy it
into this example directory as `./axisem3d`. Then launch the solver with MPI,
reading `sim/input` and writing to `sim/output`:

```bash
cd sim
rm -rf output
mpirun -np 16 axisem3d --input input --output output
```

`16` MPI ranks is plenty for the coarse 15 s mesh (use `1` on a laptop, more for
the refined meshes). On a cluster, wrap the `mpirun` command in your scheduler's
job script and load your compiler, MPI, netCDF and HDF5 modules.

Surface displacement is written in the **RTZ** frame, channel `U`, every 5th
time step (`DTx5`), in `NETCDF` format (per-rank files
`axisem3d_synthetics.nc.<rank>`) under `sim/output/stations/surface_line/`.

### 4. Figures

`make_figures.py` reads the meridian receiver line **directly** from this
example's run output `sim/output/stations/surface_line/` (the `INPUT` variable at
the top of the script; override with the `AXISEM3D_OUTPUT` environment variable,
and it prefers a higher-resolution `sim/output_1s/2s/4s` if you have one). It
maps each receiver longitude to epicentral distance Δ via `moon_line.txt` and
reconstructs the full sphere by rotational symmetry (a map point at (lat, lon)
takes the value at Δ = arccos(cos lat · cos lon)). Run it with a Python that has
`numpy`, `netCDF4`, `matplotlib`, and (for the orthographic snapshots) `cartopy`:

```bash
python make_figures.py
```

The script finds the shared `pub_style.py` automatically one directory up
(`examples/pub_style.py`), so run it from this example directory. It writes both
PDF and PNG of each figure into `figures/`:

| Figure file | Shows |
|-------------|-------|
| `moon_surface_snapshots` | 2×4 orthographic lunar-surface vertical-displacement snapshots (t ≈ 80–740 s), reconstructed from the meridian line by rotational symmetry; impact marked with a star, each panel self-scaled. |
| `moon_record_section`    | Wiggle record section: 45 normalised traces, epicentral distance (deg) vs. time after impact. |
| `moon_seismogram`        | Three-component (R, T, Z) displacement seismogram at the receiver nearest Δ ≈ 90°. |
| `moon_source_receiver_map` | Orthographic lunar globe of the source–receiver geometry (impact star + surface line, no Earth coastlines; needs cartopy). |
| `moon_overview`          | Combined figure: the first six surface snapshots (3×2) beside the record section. |

(The shipped `figures/` holds the first three as reference; running the script
regenerates all five.) If cartopy is missing, the orthographic snapshots fall
back to a plain lon/lat raster and the source–receiver map is skipped; the record
section and seismogram still build.

An optional advanced 3-D visualization (draping the wavefield over a textured
lunar globe with PyVista/Open3D) is documented in `README_gallery_pyvista_viz.md`
and uses the `conda_environment.yaml` environment plus the `stations_processing.py`
→ `png_creation_seismo.py` pipeline; it is not required to produce the figures
above.

## References

- Garcia, R. F., et al. (2019). *Lunar Seismology: An Update on Interior Structure
  Models.* Space Science Reviews — the 1-D lunar model `ISSI_MOON_M1`.
- The `.bm` model file originates from Ceri Nunn's `impact_simulations` repository:
  <https://github.com/cerinunn/impact_simulations>

## Notes

Prepared/updated by Jonathan Wolf.
