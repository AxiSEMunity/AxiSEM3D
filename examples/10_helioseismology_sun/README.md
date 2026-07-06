# Example 10: Extreme fully-acoustic simulation — solar helioseismology

An AxiSEM3D simulation of acoustic (pressure) waves propagating through the **Sun**,
modeled as a spherically symmetric, **fully fluid** body. A fluid-pressure point source
near the north pole excites the wavefield, which is recorded on a distance × azimuth grid
of 57 receivers from the pole to the antipode. Because every point in the solar model is
fluid (no shear strength, `vs = 0` everywhere), the entire mesh is acoustic — there are no
solid points. At production resolution the mesh has ~1.08 million elements and ~17.3 million
fluid points, which is what makes this an *extreme* run; it also writes a full meridional
cross-section for a whole-Sun wavefield movie. The setup is the helioseismology analogue of
a global seismic-wave simulation on Earth.

## Contents

```
10_helioseismology_sun/
├── README.md                       # this file
├── make_figures.py                 # three station figures: receiver map, seismograms, record section
├── make_cross_section.py           # interior-wavefield overview (needs the element output)
├── input/                          # all solver inputs
│   ├── solar_master.bm             # 1-D solar background model (radius, rho, vp, vs=0)
│   ├── inparam.model.yaml          # mesh + geodesy + (deactivated) 3-D cij model 'isoD'
│   ├── inparam.nr.yaml             # azimuthal resolution: CONSTANT Nr = 7
│   ├── inparam.source.yaml         # time axis + FLUID_PRESSURE source 'North_Pole'
│   ├── inparam.output.yaml         # station group + meridional cross-section element group
│   ├── inparam.advanced.yaml       # advanced/runtime options
│   └── stations_for_plotting.txt   # 57 receivers (distance × azimuth grid)
└── figures/                        # reference figures produced by the plotting scripts
```

The plotting scripts import the shared `pub_style.py` helper, which is staged one
level up at `examples/pub_style.py`; both scripts add that directory to `sys.path`
relative to themselves, so no setup is needed.

The Exodus mesh itself is **not** shipped (it is far larger than 4 MB); it is generated
from `input/solar_master.bm` with the command in the [Mesh](#mesh) section. The coarser and
full-resolution meshes are also available pre-built in the gallery (see [Large data](#large-data)).

**Key parameters (from the shipped configs):**

- **Source.** A `FLUID_PRESSURE` point source `North_Pole` at
  `latitude_longitude: [89.9, 45.0]`, `depth: 50000000.0` m (placed at depth so the source
  sits where density is moderate, avoiding the `1/rho` blow-up near the photosphere),
  amplitude `data: [1e22]`. Source-time function: `GaussianSTF`, `half_duration: 1.0` s,
  `use_derivative_integral: ERF`, `decay_factor: 1.628`, `time_shift: 0`. Setting
  `geodesy.lat_lon_north_pole_mesh: SOURCE` aligns the mesh axis with the source, so the
  wavefield is essentially axisymmetric.
- **Time axis.** `record_length: 12000.` s, `enforced_dt: NONE`, `Courant_number: 0.4`,
  `integrator: SYMPLECTIC` (energy-conserving, chosen for this long, instability-prone run).
- **Azimuthal resolution.** `inparam.nr.yaml`: `type_Nr: CONSTANT`, `constant: 7`,
  `bound_Nr_by_inplane: true`. A single near-axial pressure source needs only a few Fourier modes.
- **1-D run, no NetCDF inputs.** `inparam.model.yaml` lists one 3-D model, `isoD`
  (a `StructuredGridV3D` `cij` anisotropy model that would read `smooth_all.nc`), but it is
  set `activated: false`. The simulation is therefore effectively 1-D and needs no extra
  NetCDF model files — the Exodus mesh is the only model input.
- **Receivers.** 57 stations in `input/stations_for_plotting.txt`, laid out as a
  distance × azimuth grid: distances `d0 … d180` (latitudes 89.9999° down to −89.9999°)
  crossed with azimuths `a0, a30, a60` (longitudes 0°, 30°, 60°). Recorded as one station
  group `stationsforplotting`: `medium: FLUID`, `channels: [U]` (displacement),
  `coordinate_frame: RTZ`, `format: ASCII_STATION`, `sampling_period: 0.05` s.
- **Cross-section output.** An element group `sun_cross_section` records `U3` (radial
  displacement) on a full meridional slice (co-latitude `0 … pi`, radius `0 … 7.0e8` m) at
  `sampling_period: 50.0` s for a whole-Sun wavefield movie.

## Mesh

The shipped 1-D model is `input/solar_master.bm` — a `bm` background model giving
`radius, rho, vp, vs` for the Sun, with `vs = 0` at every depth (a purely fluid body).
The Exodus mesh referenced by `inparam.model.yaml` (`exodus_mesh: solar_master_epw4.e`) is
generated from this `.bm` file with the Salvus mesh-lite AxiSEM interface at **mesh period
500 s** and **4 elements per wavelength**:

```bash
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model solar_master.bm \
    --basic.period 500 \
    --advanced.elements_per_wavelength 4 \
    --output_filename solar_master_epw4.e
```

Run this in `input/` (so the output `solar_master_epw4.e` lands next to the configs). The
resulting mesh has ~1,080,586 elements and a solver minimum period of ~382 s.

A lighter mesh suitable for a quick test is produced from the same model at the default
element count (drop the `--advanced.elements_per_wavelength` flag), giving
`solar_new_test_500.e` (~0.27 M elements, minimum period ~320 s). If you use the lighter
mesh, set `exodus_mesh: solar_new_test_500.e` in `inparam.model.yaml`.

## Large data

This example needs **no external 3-D model data** — the only 3-D model is deactivated, so
the Exodus mesh is the sole model input. The mesh is the only large file; you can either
generate it from `solar_master.bm` with the command above, or download a pre-built copy
from the gallery:

| File | Period / EPW | Approx. size | Used by |
|------|--------------|--------------|---------|
| `solar_master_epw4.e` | 500 s, EPW 4 | ~38 MB | the shipped config (production) |
| `solar_new_test_500.e` | 500 s, default EPW | ~9.5 MB | lighter quick test |

Both are derived purely from `solar_master.bm`; no other data is required.

## Reproduce

This example assumes you have built AxiSEM3D and have the `axisem3d` binary on your
`PATH` (or copied it into this directory as `./axisem3d`). The solver CLI flags are
`--input` / `--output` (not `-i` / `-o`).

### 1. Mesh

Generate the Exodus mesh from `input/solar_master.bm` with the `salvus_mesh_lite`
command in the [Mesh](#mesh) section, placing the result in `input/`. The shipped config
expects `input/solar_master_epw4.e` (EPW 4); for a quick test use the lighter default-EPW
mesh and set `exodus_mesh: solar_new_test_500.e` in `input/inparam.model.yaml`.

### 2. Large data

None beyond the mesh — the only 3-D model is deactivated. The mesh can be generated
locally (step 1) or downloaded pre-built from the gallery (see [Large data](#large-data)).

### 3. Run

With the mesh in `input/`, launch the solver with MPI from this example directory:

```bash
mpirun -np <N> axisem3d --input input --output output
```

On a cluster, wrap the `mpirun` command in your scheduler's job script and load your
compiler, MPI, netCDF and HDF5 modules.

Choose `<N>` for your hardware: the lighter default-EPW mesh (~0.27 M elements) runs
comfortably on a few dozen ranks (e.g. `-np 24`); the production EPW4 mesh (~1.08 M
elements) wants hundreds (e.g. `-np 240`). With the EPW4 mesh and `Courant_number: 0.4`,
the solver picks
Δt ≈ 0.0411 s and takes **292,104** SYMPLECTIC time steps for the 12,000 s record.
Output is written under `output/`:

- `output/stations/stationsforplotting/<network>.<name>.ascii` — one ASCII file per receiver
  (57 total; the naming is `<network>.<name>`, e.g. `a0.d0.ascii`). Each file has three
  columns = the `RTZ` displacement components.
- `output/stations/stationsforplotting/data_time.ascii` — the shared time axis (seconds).
- `output/elements/sun_cross_section/axisem3d_synthetics.nc.rank*` — the meridional
  cross-section element output (per-rank NetCDF), used for the interior-wavefield figure.
- `output/develop/preloop_diagnosis.log` — solver diagnostics.

### 4. Figures

Two plotting scripts build the reference figures into `figures/` (each as `.pdf` + `.png`).
Both read this example's own `output/` relative to themselves and import the shared
`pub_style.py` from the examples root, so they run from anywhere with no path edits.

```bash
python make_figures.py        # needs numpy, pyyaml, matplotlib (no cartopy)
python make_cross_section.py  # additionally needs pandas and netCDF4
```

**`make_figures.py`** reads the ASCII station output
(`output/stations/stationsforplotting/`) and the source/receiver lists in `input/`, and
writes three figures:

- **`sun_receiver_map`** — a polar map centered on the pole source: each receiver is a
  triangle at (azimuth = longitude, radius = great-circle angular distance from the source);
  the source is a star at the center.
- **`sun_seismograms`** — the three-component (R, T, Z) displacement seismogram at a
  mid-distance receiver, with the epicentral distance in the title.
- **`sun_record_section`** — a vertical-displacement record section along the `a0` azimuth
  line: each trace is amplitude-normalized and offset by its epicentral distance, showing the
  acoustic wavefront move out from pole toward antipode.

**`make_cross_section.py`** reads the meridional element wavefield from
`output/elements/sun_cross_section/` (raw per-rank NetCDF) and the `a0` station line, and
writes one combined overview figure:

- **`sun_overview`** — left: the first six interior-wavefield snapshots (radial displacement
  `U3`, mirrored about the rotation axis, each normalized to its own peak), showing the
  acoustic wavefront sweep through the solar interior; right: the `a0`-line vertical-
  displacement record section. This figure requires the element cross-section group, so it
  is only produced after running the full simulation (step 3).

## Notes

- **All-fluid mesh.** Every point in the solar model is a fluid point (`vs = 0`), so
  `SolidPoint$Mass1D` never appears in the run log and station groups recording in this mesh
  must use `medium: FLUID`.
- **No NetCDF model inputs.** The only 3-D model (`isoD`, a `cij` anisotropy model) is
  `activated: false`, so `smooth_all.nc` is never read; the run is effectively 1-D.
- **Stability.** The combination of `SYMPLECTIC` integration and `Courant_number: 0.4`
  (rather than the default ~0.6) keeps this long acoustic run stable; placing the source at
  depth avoids the `1/rho` amplitude blow-up near the photosphere.

Prepared/updated by Jonathan Wolf.
