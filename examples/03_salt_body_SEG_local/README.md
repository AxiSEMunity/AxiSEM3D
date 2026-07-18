# Example 03: Local SEG/EAGE C3 salt body (5 Hz)

**Author:** Jonathan Wolf (wolf@ucsc.edu)

A local, high-frequency exploration-seismology simulation on a 2D-axisymmetric
Cartesian mesh. An explosive source just below an ocean floor radiates through
the **SEG/EAGE C3 salt-body model** (the canonical 3D salt benchmark). A thin
water layer (fluid) sits on top of a solid sediment/salt section, so the
recording surface is the **ocean floor**. A 1D reference run (`input1D`) is
compared with the full 3D run (`input3D`), in which the volumetric velocity model
and the undulating salt-body interface are read from NetCDF data files. The
ocean-floor waveforms and the in-plane wavefield cross-sections show the imprint
of the salt body.

The simulation runs at its native **0.2 s dominant period (~5 Hz)** — it is not
rerun at the longer periods used by the global examples.

## Contents
- `SEG_C3.bm` — polynomial background (1D reference) model used to build the mesh.
- `input1D/` — input parameter files (`inparam.model/nr/source/output/advanced.yaml`)
  for the 1D run, plus the mesh.
- `input3D/` — input parameter files for the 3D run, plus the mesh.
  `inparam.model.yaml` adds the volumetric model (`SEG_C3_data/SEG_C3_SOLID.nc`)
  and the geometric salt-body undulation (`SEG_C3_data/SEG_C3_UNDULATION.nc`),
  both referenced by relative name.
- `local_mesh__SEG_salt__5Hz.e` — the shipped coarse Exodus mesh (present in both
  input folders).
- `make_figures.py` — reproduces the figures in `figures/`.
- `post_processing.ipynb` — notebook that reads the element output and plots
  waveforms / writes VTK slices for animation.
- `data/salt_model.npz` — a compact (~0.4 MB) extract of the input velocity model,
  used only for the optional salt-model figure (see **Figures** below).
- `figures/` — reference figures.

## Solver binary

The run commands below assume you have already built AxiSEM3D and that the solver
is available as `axisem3d` — either on your `PATH`, or copied into this example
directory as `./axisem3d`. Nothing here hardcodes a build path; adjust the MPI
rank count (`-np`) to suit your machine. On a cluster, wrap the `mpirun` command
in your scheduler's job script and load your compiler, MPI, netCDF and HDF5
modules.

## Reproduce

### 1. Mesh
The shipped mesh `local_mesh__SEG_salt__5Hz.e` is a 5 Hz (0.2 s period) Cartesian
AxiSEM mesh built with `salvus_mesh_lite` from the polynomial background model
`SEG_C3.bm` (two regions: a 200 m water layer with VS = 0 over solid sediment,
VP ≈ 1.5–1.6 km/s). The domain is 5 km wide; the top of the model is r = 6371 km,
the ocean floor at r = 6370.8 km. To (re)generate it:
```bash
python -m salvus_mesh_lite.interface AxiSEMCartesian \
    --basic.model SEG_C3.bm --basic.period .2 \
    --cartesian2Daxisem.x 5. --cartesian2Daxisem.min_z 6367.0 \
    --attenuation.frequencies 0.01 10. \
    --output_filename local_mesh__SEG_salt__5Hz.e
```

### 2. Large data (3D run only)
The 3D run needs the SEG/EAGE C3 model data referenced by
`input3D/inparam.model.yaml`:
- `SEG_C3_data/SEG_C3_SOLID.nc` — volumetric VP/VS/RHO model (~1 GB), applied as a
  `StructuredGridV3D` model on an `XY_CARTESIAN` grid.
- `SEG_C3_data/SEG_C3_UNDULATION.nc` — geometric undulation (0–1000 m) of the
  200 m interface (~1.8 MB), applied as a `StructuredGridG3D` model.

These exceed the repository size limit and are **not** shipped here. Obtain the
archive `SEG_C3_data.tar.bz2` from the example gallery / Zenodo release, place it
in `input3D/`, and extract it there so the model files sit at the relative paths
`input3D/inparam.model.yaml` expects:
```bash
tar -xjf SEG_C3_data.tar.bz2 -C input3D
# result: input3D/SEG_C3_data/SEG_C3_SOLID.nc and
#         input3D/SEG_C3_data/SEG_C3_UNDULATION.nc
```
The 1D run needs no extra data.

### 3. Run
The source is an explosive `MOMENT_TENSOR` at depth 1000 m, with a Gaussian
source-time function (half-duration 0.2 s), record length 10 s, Courant 0.5,
NEWMARK, attenuation off. Output is element-wise NetCDF (two element groups; see
**Output** below). On the shipped coarse mesh, 4 MPI ranks is plenty. Run each
case directly against its input folder, writing the results to a matching
`simu*/output/` directory:

```bash
# 1D reference (~0.3 min on 4 cores)
OMP_NUM_THREADS=1 mpirun -np 4 axisem3d --input input1D --output simu1D/output

# full 3D salt body (~10 min on 4 cores; needs the SEG_C3 data from step 2)
OMP_NUM_THREADS=1 mpirun -np 4 axisem3d --input input3D --output simu3D/output
```

`make_figures.py` reads `simu3D/output/` by default, so keeping the output paths
above lets the figure step work out of the box. (You may write elsewhere and
point the figure script at it with the `EX03_OUTPUT` variable — see **Figures**.)

### 4. Figures
```bash
python make_figures.py
```
`make_figures.py` reads the 3D element output under `simu3D/output/elements/`
(override the location with the `EX03_OUTPUT` environment variable) and writes
both PDF and PNG to `figures/`:

1. **`ex03_ocean_floor_waveforms`** — 3 stacked panels (Radial / Transverse /
   Vertical, in metres, shared y-scale) comparing three ocean-floor locations:
   (2, 3) km, (5, 0) km, (0, 5) km. The time series at each (x, y) is
   reconstructed from the stored azimuthal Fourier coefficients. The Gaussian STF
   is acausal, so the time axis starts at t = 0.
2. **`ex03_cross_section_snapshots`** — an 8-panel (2×4) time-lapse of the
   vertical displacement $U_Z$ on the φ = 0 in-plane slice, drawn with
   `tricontourf` over the triangulated GLL points (shared colour scale, one
   colorbar). The y-axis is depth below the top of the imaged model (m, increasing
   downward); snapshots span ≈8–85 % of the record.
3. **`ex03_salt_model_velocities`** *(optional)* — the **input** velocity model of
   the SEG/EAGE C3 salt body (not a wavefield): a 2×2 panel showing a vertical
   cross-section and a horizontal map slice for each of $V_P$ (top, *viridis*) and
   $V_S$ (bottom, *magma*). The salt body is the high-velocity blob
   ($V_P$ ≈ 4.5 km/s, $V_S$ = 3 km/s) inside the slow layered sediment.

   This figure does **not** come from running the example — it visualises the
   input model. The script reads `data/salt_model.npz`, a compact (~0.4 MB)
   extract of two slices (one vertical, one horizontal) of the full ~1 GB
   `SEG_C3_SOLID.nc` input model. That `.npz` is shipped here so the figure works
   out of the box; if you delete it, the script skips Fig 3. To regenerate it from
   the source `.nc`, sample `VP`/`VS` along the salt-richest vertical y-slice and
   the widest-extent horizontal depth-slice and save the arrays
   (`x, y, depth, ybest, dbest, VP_xz, VS_xz, VP_xy, VS_xy`) with
   `numpy.savez_compressed`.

Shared styling (colorblind-safe RTZ palette, serif fonts) comes from
`examples/pub_style.py`, one directory up; `make_figures.py` adds that directory
to `sys.path` relative to itself, so no installation or `PYTHONPATH` tweak is
needed.

## Output
Element-wise NetCDF output (no station groups), two element groups defined in
`inparam.output.yaml`, both displacement `U` in the RTZ frame, written per MPI
rank as `axisem3d_synthetics.nc.rankN`:

- **`Fourier_coefficients_ocean_floor`** — the vertical edge at the ocean floor
  (`edge_position: 6370.8e3`), all GLL points, **full Fourier coefficients**
  (`phi_list: []`), sampling period 0.01 s. Used for the waveform figure.
- **`orthogonal_azimuthal_slices`** — whole elements (`edge_dimension: BOTH`, GLL
  subset `[0, 2, 4]`), four azimuthal slices `phi_list = [0, π/2, π, 3π/2]`,
  sampling period 0.05 s. Used for the cross-section time-lapse.

## Post-processing notebook
`post_processing.ipynb` reads the element output under `simu3D/output/elements/`
and shows how to plot the ocean-floor waveforms and how to export the in-plane
slices to VTK for animation.
