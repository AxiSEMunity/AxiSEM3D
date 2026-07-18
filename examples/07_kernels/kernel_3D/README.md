# Example 07: 3-D high-frequency sensitivity kernels — kernel_3D (S362ANI + Moho topography)

P-wave traveltime (banana-doughnut) and Moho interface sensitivity kernels computed
from an AxiSEM3D **3-D** forward simulation — S362ANI volumetric heterogeneity plus a
Gaussian Moho-topography perturbation — using the external
[`axikernels`](https://github.com/Adrian-Mag/AxiSEM3D_Kernels) package. It is the 3-D
counterpart of [`../kernel_1D/`](../kernel_1D/): the same vp and Moho (P400P) measurements,
but on a laterally heterogeneous Earth instead of plain PREM. The forward run dumps the
full elastic wavefield over a pole-to-pole element region; `axikernels` then builds the
adjoint source, runs the backward simulation, and evaluates the kernels. The Moho interface
kernel `K_SS` spatially correlates with the imposed Gaussian topography (overlaid as contours
in the figure). Kernel computation is **not** built into AxiSEM3D core.

## Contents

```
kernel_3D/
  README.md                  -- this file
  readme.txt                 -- upstream axikernels companion notes (background, API)
  compute_kernels_3D.py      -- axikernels driver: adjoint source -> backward run -> vp + Moho kernels
  convert_s362ani_to_radius.py -- regenerates input_forward/S362ANI_radius.nc (see Large data)
  create_moho_topography.py  -- regenerates input_forward/moho_topography.nc (Gaussian Moho perturbation)
  make_figures.py            -- optional re-styling of the kernel HDF5 into publication figures
  input_forward/             -- AxiSEM3D inputs for the forward run
    inparam.model.yaml         -- 1-D PREM mesh + 3-D model list (Moho topography, then S362ANI)
    inparam.source.yaml        -- isotropic moment-tensor (explosion) source, GaussianSTF 25 s, record_length 600 s
    inparam.nr.yaml            -- azimuthal resolution: CONSTANT Nr = 50 (resolves the 3-D structure)
    inparam.output.yaml        -- element output, channels [U, G, S], RTZ
    inparam.advanced.yaml      -- advanced solver options
    prem_iso_elastic.bm        -- PREM 1-D reference / background model card (Moho at radius 6 346 600 m)
    prem_iso_elastic_50s.e     -- 50 s Exodus mesh (used by the run)
    moho_topography.nc         -- Gaussian Moho topography (vars latitude/longitude/undulation_MOHO)
    STA_10DEG_GRID.txt         -- 10-degree station grid
    (S362ANI_radius.nc)        -- volumetric S362ANI model; NOT shipped, see Large data
  figures/                   -- reference figures (kernel3D_vp_slice.{pdf,png}, kernel3D_moho.{pdf,png})
```

**3-D model** (`inparam.model.yaml`, `list_of_3D_models`, applied in this order):

1. `MOHO_TOPOGRAPHY` — a Gaussian Moho-depth perturbation about the PREM Moho radius
   (`interface: 6346600.` m, depth ≈ 24.4 km), read from `moho_topography.nc` as a
   `StructuredGridG3D` interface model (the mesh geometry is undulated first).
2. `EMC_S362ANI` — S362ANI volumetric VP/VS heterogeneity as a `StructuredGridV3D` model
   from `S362ANI_radius.nc`, applied as a percentage perturbation to the 1-D reference
   (`reference_kind: REF1D`, `factor: 5.0` — amplitude exaggeration to make the kernels
   visually obvious).

The source matches `kernel_1D` (isotropic moment tensor `[1e21, 1e21, 1e21, 0, 0, 0]` N·m,
`latitude_longitude: [0, 0]`, depth 20 km, `GaussianSTF` 25 s, `record_length: 600.` s,
`Courant_number: 0.7`, `NEWMARK`). The key difference is `inparam.nr.yaml`:
`type_Nr: CONSTANT`, `constant: 50` — a 3-D model (S362ANI, spherical-harmonic degree 18)
at 50 s needs Nr ≥ 50 to resolve the lateral heterogeneity (the 1-D example uses Nr = 5).

## Reproduce

You need the `axisem3d` solver: build AxiSEM3D and put the `axisem3d` binary on your
`PATH`, or copy it into this directory as `./axisem3d` (`cp /path/to/axisem3d .`). The
kernel step also needs `axikernels` installed (see `readme.txt`). On a cluster, wrap the
`mpirun` command below in your scheduler's job script and load your compiler, MPI, netCDF
and HDF5 modules. The steps below run in order.

1. **Mesh** — the 50 s Exodus mesh `input_forward/prem_iso_elastic_50s.e` ships; regenerate
   it from `prem_iso_elastic.bm` only if you want a different period (see *Mesh*).
2. **Large data** — the volumetric `input_forward/S362ANI_radius.nc` (~9.4 MB) is **not
   shipped** (exceeds the size limit). Regenerate it with `python convert_s362ani_to_radius.py`
   or download it from the gallery, and regenerate the Gaussian Moho file with
   `python create_moho_topography.py` if needed (see *Large data*). Both NetCDF files must
   be present in `input_forward/` before the forward run.
3. **Run the forward simulation.** The kernel evaluation reads the FORWARD run output as
   its input, so the forward solver runs first. The `axikernels` driver expects a single
   simulation directory containing `input/`, `output/`, and the `axisem3d` binary (which it
   reuses for the backward run), so assemble `simu_forward/` and run the solver inside it:

   ```bash
   mkdir -p simu_forward
   cp -r input_forward simu_forward/input          # includes the NetCDF model files
   cp /path/to/axisem3d simu_forward/axisem3d      # or your built binary
   mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output
   ```

   The pole-to-pole `[U, G, S]` element dump lands in `simu_forward/output/elements/`. Then
   build the adjoint source, run the backward simulation, and evaluate the kernels:

   ```bash
   python compute_kernels_3D.py --cores 8 --forward simu_forward \
       --output kernel_output --topography input_forward/moho_topography.nc
   ```

   This **creates and runs the backward simulation automatically** from the forward output
   (a sibling directory `backward_simu_forward/`, launched internally as
   `mpirun -np <cores> axisem3d`), writing `kernel_output/{vp_kernel,moho_kd}.{h5,png}`.
   `--forward` / `--output` / `--topography` default to `simu_forward` / `kernel_output` /
   `input_forward/moho_topography.nc`, so bare `python compute_kernels_3D.py` does the same
   thing. No manual backward run is needed.
4. **Figures** — re-style the kernel HDF5 into the publication figures:

   ```bash
   python make_figures.py
   ```

   Reads `kernel_output/` (the `INPUT` path variable at the top of the script; override
   with the `KERNEL_OUTPUT` environment variable) and overlays the imposed Moho topography
   from `input_forward/moho_topography.nc`, writing to `figures/`:

   | Figure | Shows |
   |--------|-------|
   | `figures/kernel3D_vp_slice.{pdf,png}` | $v_P$ traveltime (banana-doughnut) kernel on the great-circle slice, 3-D model |
   | `figures/kernel3D_moho.{pdf,png}`     | Moho interface kernel `K_SS` with imposed-topography contours overlaid |

   `make_figures.py` imports the shared `pub_style.py` from the examples root (resolved
   relative to the script) and uses the matplotlib `Agg` backend, so it runs headless.

## Mesh

The run uses the **50 s** coarse Exodus mesh `input_forward/prem_iso_elastic_50s.e`
(referenced by `inparam.model.yaml` as `exodus_mesh: prem_iso_elastic_50s.e`), which
ships in `input_forward/`. Regenerate it from the 1-D model card with the AxiSEM3D mesher:

```bash
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model prem_iso_elastic.bm --basic.period 50.0 \
    --output_filename prem_iso_elastic_50s.e
```

## Large data

The volumetric model file **`input_forward/S362ANI_radius.nc`** (~9.4 MB, S362ANI on a
radius grid) exceeds the 4 MB limit and is **not shipped here**. Obtain it one of two ways:

- **Regenerate it** with the included script, which reads `S362ANI_percent.nc` from the
  axikernels example data (`AxiSEM3D_Kernels/examples/data/3D_KERNEL_EXAMPLE_30s/input/`)
  and converts depth→radius:

  ```bash
  python convert_s362ani_to_radius.py        # -> input_forward/S362ANI_radius.nc
  ```

- **Or download it** from the example gallery (see the gallery manifest accompanying this
  example) and place it in `input_forward/`.

The Gaussian Moho file `input_forward/moho_topography.nc` (~0.5 MB) ships, and can be
regenerated with `python create_moho_topography.py`. The 3-D forward and backward element
wavefield dumps are very large (multiple GB each) and are **generated** by the run, not
shipped.

## Run details

The validated kernel parameters (also the script defaults) are `--receiver 0 40`,
`--window 425 475`, `--tau 2`, `--channel UZ`, `--resolution 200`; the forward and
backward runs each use 8 MPI ranks (`-np 8` / `--cores 8`).

Outputs in `kernel_output/`: `vp_kernel.{h5,png}` and `moho_kd.{h5,png}`. The `.png` files
are produced directly by `compute_kernels_3D.py` and are the minimum reproducible figures;
`make_figures.py` (step 4) is the optional re-styling into the publication figures under
`figures/` (and overlays the Moho topography from `input_forward/moho_topography.nc`). To
start fresh, delete the generated `simu_forward/`, `backward_simu_forward/`, and
`kernel_output/` directories (all git-ignored).

## Notes

- Companion example: kernel computation is done by the external `axikernels` package on
  the AxiSEM3D wavefields, not in AxiSEM3D core. See `readme.txt` for the package and its
  environment.
- **Model order matters.** `MOHO_TOPOGRAPHY` is listed before `EMC_S362ANI` so the mesh
  geometry is undulated first, then the volumetric model is applied.
- **Nr is 50, not 5.** The 3-D S362ANI model needs `constant: 50` to resolve the degree-18
  lateral structure at 50 s; lowering it under-resolves the heterogeneity.
- **Moho radius.** The Moho is the PREM solid-solid discontinuity at r = 6 346 600 m
  (depth ≈ 24.4 km), set in `inparam.model.yaml` (`undulation_range.interface`) and
  confirmed from `prem_iso_elastic.bm`. AxiSEM3D stores 3-D element output in reference
  (undeformed-mesh) coordinates, so the same Moho radius applies to both the 1-D and 3-D runs.

## References

- **AxiSEM3D** — Leng, Nissen-Meyer, van Driel et al., *Geophys. J. Int.* (2019);
  <https://github.com/AxiSEMunity/AxiSEM3D>.
- **axikernels** — <https://github.com/Adrian-Mag/AxiSEM3D_Kernels> (Adrian Mag).
- **S362ANI** — Kustowski, Ekström & Dziewoński (2008), *J. Geophys. Res.* (via IRIS EMC).
- **PREM** — Dziewoński & Anderson (1981), *Phys. Earth Planet. Inter.* 25, 297–356.
- Companion 1-D example: [`../kernel_1D/`](../kernel_1D/).

Prepared/updated by Jonathan Wolf.
