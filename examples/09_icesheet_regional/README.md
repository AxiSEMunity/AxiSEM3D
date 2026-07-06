# Example 09: Complex multi-media propagation — the Enceladus ice shell

A regional-scale AxiSEM3D simulation of seismic wave propagation through the
ice shell of Saturn's moon Enceladus. The model couples multiple media — a
solid ice crust, a subsurface liquid-water ocean, and an underlying silicate
upper layer — and includes undulating internal and surface topography, so it
exercises AxiSEM3D's handling of solid/fluid coupling, geometric (boundary)
3D models, and volumetric 3D models in a single run. The 1D background model
and the volumetric ice-shell properties follow
[Dapré et al., 2024](https://doi.org/10.1029/2024JE008644), and the basal
(ice–ocean) topography follows
[Cadek et al., 2016](https://doi.org/10.1002/2016GL068634). The source is a
contractional moment tensor at the sub-source point (0°, 0°), 2 km below the
surface, recorded by a regional network of 81 surface stations. Because this is
a genuinely 3D model with a moment-tensor source, the transverse component
carries energy (SH / Love-type arrivals), unlike a purely axisymmetric run.

## Contents
- `input/inparam.model.yaml` — 1D mesh, geodesy, absorbing boundaries,
  attenuation, and the two 3D models: volumetric ice-shell properties
  (`StructuredGridV3D`) and the Cadek basal topography (`StructuredGridG3D`).
- `input/inparam.source.yaml` — time axis (record length 1000 s, Courant 0.6,
  NEWMARK) and the moment-tensor source with a Gaussian source-time function.
- `input/inparam.nr.yaml` — azimuthal (Fourier) resolution: `CONSTANT`, Nr = 5.
- `input/inparam.output.yaml` — station group `JW_stations`, RTZ displacement
  (`channels: [U]`), `ASCII_STATION` format.
- `input/inparam.advanced.yaml` — advanced/runtime settings.
- `input/Enceladus.bm` — 1D radial background model (used to build the mesh).
- `input/enceladus_2s_chunk.e` — shipped coarse Exodus mesh (2 s period; see
  **Mesh**).
- `input/Cadek_basal_topography.nc` — basal ice–ocean topography grid.
- `input/3Dstations.txt` — station network (name, network, latitude, longitude,
  depth); 81 surface stations on a lat/lon grid.
- `make_figures.py` — post-processing script that produces the figures.
- `figures/` — reference figures (record section, station map, seismograms,
  combined overview).

## Output
Station output is written in `ASCII_STATION` format (one file per station, all
channels in columns) in the **RTZ** frame, so the columns are
`U1 = R, U2 = T, U3 = Z` — displacement in metres. The time vector is in
`data_time.ascii` (seconds after source origin). Files land in
`output/stations/JW_stations/`.

## Mesh
The shipped mesh `input/enceladus_2s_chunk.e` is a **2 s** (minimum period)
axisymmetric chunk covering a 60° colatitude cap. Regenerate or refine it from
the 1D background model `Enceladus.bm` with the salvus mesher (any machine with
`salvus_mesh_lite` on its Python path):

```bash
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model Enceladus.bm \
    --basic.period 2 \
    --spherical.min_radius 150 \
    --chunk2D.max_colatitude 60 \
    --output_file enceladus_2s_chunk.e
```

Change `--basic.period` to mesh at a different minimum period (a smaller period
gives a finer, larger mesh) and update `exodus_mesh:` in `inparam.model.yaml` to
the new filename.

## Reproduce

This example assumes you have built AxiSEM3D and either have the `axisem3d`
binary on your `PATH` or have copied it into this directory as `./axisem3d`.

1. **Mesh.** The 2 s mesh `input/enceladus_2s_chunk.e` ships with the example —
   no action needed to reproduce the shipped run. To refine, regenerate it as in
   **Mesh** above.

2. **Large data.** The volumetric ice-shell model
   `volumetric_ice_shell_properties.nc` (`nc_data_file` in `inparam.model.yaml`,
   ~44 MB) is too large to ship in the repository and is **not** included here.
   Obtain it from the example gallery and place it in `input/` before running.
   It follows [Dapré et al., 2024](https://doi.org/10.1029/2024JE008644).

3. **Run.** From this example directory, launch the solver under MPI, pointing
   it at `input/` and `output/`:

   ```bash
   mpirun -np 6 axisem3d --input input --output output
   ```

   `-np 6` suits the shipped coarse 2 s mesh; scale the rank count to your
   machine and mesh resolution. Use `./axisem3d` if you copied the binary into
   this directory rather than putting it on your `PATH`. On a cluster, wrap the
   `mpirun` command in your scheduler's job script and load your compiler, MPI,
   netCDF and HDF5 modules.

   The station seismograms are written to `output/stations/JW_stations/`.

4. **Figures.** Generate the figures with:

   ```bash
   python make_figures.py
   ```

   `make_figures.py` reads the run output from `output/stations/JW_stations/`
   (set by the `INPUT` variable at the top of the script) together with the
   station-info file `input/3Dstations.txt`, and imports the shared plotting
   helper `pub_style.py` from the top-level `examples/` directory (one level
   above this example). It needs `numpy` and `matplotlib`; the combined overview
   figure additionally uses `scipy` (and `cartopy` for the map projection, if
   installed — it falls back to a plain lat/lon grid otherwise). The figures are
   written to `figures/` as PNG and PDF:

   - **`ex09_enceladus_record_section`** — vertical-component (Z) record section,
     traces **sorted by epicentral distance** and offset by distance so the
     moveout across the receiver grid is visible.
   - **`ex09_enceladus_station_map`** — lat/lon scatter of the full station
     network (blue triangles) with the **source marked as a star**, equal aspect.
   - **`ex09_enceladus_seismograms`** — three-component (R/T/Z) seismogram at the
     station where the **transverse arrives earliest relative to the radial**
     (per-component onset detection, so a weak-but-early SH arrival is caught);
     R/T/Z share a single y-scale (metres) so relative amplitudes are honest, and
     the epicentral distance Δ is printed in the title.
   - **`ex09_enceladus_overview`** — combined overview: six snapshots of the
     surface vertical wavefield, interpolated across the real station grid (each
     station normalised to its own peak so the travelling wavefront is visible),
     beside the vertical-component record section.

   Colours are colorblind-safe (Okabe-Ito, from `pub_style.py`).

## Notes
Prepared by Jonathan Wolf.
