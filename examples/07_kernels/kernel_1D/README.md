# Example 07: 3-D high-frequency sensitivity kernels — kernel_1D (PREM)

P-wave (vp) cross-correlation traveltime **banana-doughnut** sensitivity kernel for
the **P400P** phase, plus a **Moho interface kernel** (`K_SS`), computed from a 1-D
PREM AxiSEM3D forward run. This is a *companion* example: AxiSEM3D is not a kernel
code — it supplies the forward (and adjoint/backward) wavefields, and the external
Python package [`axikernels`](https://github.com/Adrian-Mag/AxiSEM3D_Kernels) builds
the adjoint source, launches the backward simulation, and evaluates the kernels.
The vp kernel sliced on the source–receiver great circle shows the classic
banana-doughnut pattern (zero sensitivity on the geometric ray, peak in the first
Fresnel zone); the Moho interface kernel reuses the same wavefields, so it costs no
extra simulation. The 3-D counterpart on a heterogeneous Earth model is in
[`../kernel_3D/`](../kernel_3D/).

## Contents

```
kernel_1D/
  README.md             -- this file
  readme.txt            -- upstream axikernels companion notes (scientific background, API)
  compute_kernels.py    -- axikernels driver: builds adjoint source, runs backward sim, evaluates kernels
  make_figures.py       -- optional re-styling of the kernel HDF5 into publication figures
  input_forward/        -- AxiSEM3D inputs for the forward run
    inparam.model.yaml    -- 1-D PREM iso, mesh prem_iso_elastic_50s.e, CG4 attenuation, [RIGHT,BOTTOM] absorbing
    inparam.source.yaml   -- isotropic moment-tensor (explosion) source, GaussianSTF 25 s, record_length 600 s
    inparam.nr.yaml       -- azimuthal resolution: CONSTANT Nr = 5 (single axial source)
    inparam.output.yaml   -- element output, channels [U, G, S] (needed to reconstruct strain)
    inparam.advanced.yaml -- advanced solver options
    prem_iso_elastic.bm   -- PREM 1-D model card (defines the Moho at radius 6 346 600 m)
    prem_iso_elastic_50s.e-- 50 s Exodus mesh (used by the run)
    prem_iso_elastic_40s.e-- 40 s Exodus mesh (alternative)
    prem_iso_elastic_30s.e-- 30 s Exodus mesh (alternative)
    STA_10DEG_GRID.txt    -- 10-degree station grid
  input_backward/       -- placeholder only; axikernels generates the backward input itself
  figures/              -- reference figures (kernel1D_vp_slice.{pdf,png}, kernel1D_moho.{pdf,png})
```

The source is a single on-axis **isotropic moment tensor** (`data: [1e21, 1e21, 1e21, 0, 0, 0]`
N·m, `latitude_longitude: [0, 0]`, depth 20 km), `GaussianSTF` with `half_duration: 25.` s.
An on-axis isotropic source radiates an axisymmetric (m = 0) field, so `Nr` can be small
(`CONSTANT = 5`). Time axis: `record_length: 600.` s, `Courant_number: 0.7`, `NEWMARK`.

## Reproduce

You need the `axisem3d` solver: build AxiSEM3D and put the `axisem3d` binary on your
`PATH`, or copy it into this directory as `./axisem3d` (`cp /path/to/axisem3d .`). The
kernel step also needs `axikernels` installed (see `readme.txt`). On a cluster, wrap the
`mpirun` command below in your scheduler's job script and load your compiler, MPI, netCDF
and HDF5 modules. The steps below run in order.

1. **Mesh** — the 50 s Exodus mesh `input_forward/prem_iso_elastic_50s.e` ships with the
   example; regenerate it from `prem_iso_elastic.bm` only if you want a different period
   (see *Mesh*).
2. **Large data** — none needed for the 1-D run; only the shipped mesh and `.bm` card are
   used.
3. **Run the forward simulation.** The kernel evaluation reads the FORWARD run output as
   its input, so the forward solver runs first. The `axikernels` driver expects a single
   simulation directory containing `input/`, `output/`, and the `axisem3d` binary (which it
   reuses for the backward run), so assemble `simu_forward/` and run the solver inside it:

   ```bash
   mkdir -p simu_forward
   cp -r input_forward simu_forward/input
   cp /path/to/axisem3d simu_forward/axisem3d        # or your built binary
   mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output
   ```

   The forward element wavefield lands in `simu_forward/output/elements/`. Then build the
   adjoint source, run the backward simulation, and evaluate the kernels with the axikernels
   driver:

   ```bash
   python compute_kernels.py --cores 8 --forward simu_forward --output kernel_output
   ```

   This **creates and runs the backward simulation automatically** from the forward output
   (a sibling directory `backward_simu_forward/`, launched internally as
   `mpirun -np <cores> axisem3d`), writing `kernel_output/{vp_kernel,moho_kd}.{h5,png}`.
   `--forward` / `--output` default to `simu_forward` / `kernel_output`, so bare
   `python compute_kernels.py` does the same thing. No manual backward run is needed, and
   `input_backward/` is never populated by the user.
4. **Figures** — re-style the kernel HDF5 into the publication figures:

   ```bash
   python make_figures.py
   ```

   Reads `kernel_output/` (the `INPUT` path variable at the top of the script; override
   with the `KERNEL_OUTPUT` environment variable) and writes to `figures/`:

   | Figure | Shows |
   |--------|-------|
   | `figures/kernel1D_vp_slice.{pdf,png}` | $v_P$ traveltime (banana-doughnut) kernel on the source–receiver great-circle slice |
   | `figures/kernel1D_moho.{pdf,png}`     | Moho interface kernel `K_SS` map |

   `make_figures.py` imports the shared `pub_style.py` from the examples root (resolved
   relative to the script) and uses the matplotlib `Agg` backend, so it runs headless.

## Mesh

The run uses the **50 s** coarse Exodus mesh `input_forward/prem_iso_elastic_50s.e`
(referenced by `inparam.model.yaml` as `exodus_mesh: prem_iso_elastic_50s.e`), which
ships in `input_forward/`. The 40 s and 30 s meshes ship as alternatives. The meshes
are generated from the 1-D model card `input_forward/prem_iso_elastic.bm` with the
AxiSEM3D mesher (`salvus_mesher`/`AxiSEM3D-Mesher`); for the 50 s mesh:

```bash
python -m salvus_mesh_lite.interface AxiSEM \
    --basic.model prem_iso_elastic.bm --basic.period 50.0 \
    --output_filename prem_iso_elastic_50s.e
```

(Substitute `--basic.period 40.0` / `30.0` to regenerate the alternative meshes.)

## Large data

None required: the 1-D run needs only the shipped mesh and `.bm` model card. The
forward and backward element wavefield dumps that the kernel evaluation consumes are
large (tens of GB combined) and are **generated** by the steps above, not shipped.

## Run details

The validated kernel parameters (also the script defaults) are `--receiver 0 40`,
`--window 425 475`, `--tau 2`, `--channel UZ`, `--resolution 200`; the forward and
backward runs each use 8 MPI ranks (`-np 8` / `--cores 8`).

Outputs in `kernel_output/`:

| File | Contents |
|------|----------|
| `vp_kernel.h5`  | vp kernel data (reloadable via axikernels) |
| `vp_kernel.png` | quick-look vp banana-doughnut slice |
| `moho_kd.h5`    | Moho interface kernel `K_SS` (+ normal/`K_dn` and material-contrast/`K_dv` parts) on a lat/lon grid; attr `moho_radius_m = 6346600.0` |
| `moho_kd.png`   | quick-look `K_SS` map |

The quick-look figures `kernel_output/vp_kernel.png` and `kernel_output/moho_kd.png`
are produced directly by `compute_kernels.py` and are the minimum reproducible figures.
`make_figures.py` (step 4) is the optional re-styling into the publication figures under
`figures/`. To start fresh, delete the generated `simu_forward/`,
`backward_simu_forward/`, and `kernel_output/` directories (all git-ignored).

## Notes

- This is a companion example: kernel computation is done by the external `axikernels`
  package on the AxiSEM3D wavefields, not in AxiSEM3D core. See `readme.txt` for the
  package, its environment, and the scientific background.
- The Moho is the PREM solid-solid discontinuity at radius **6 346 600 m** (depth 24.4 km),
  confirmed from `prem_iso_elastic.bm`.

## References

- **PREM** — Dziewoński & Anderson (1981), *Phys. Earth Planet. Inter.* 25, 297–356.
- **AxiSEM3D** — Leng, Nissen-Meyer, van Driel et al., *Geophys. J. Int.* (2019);
  <https://github.com/AxiSEMunity/AxiSEM3D>.
- **axikernels** — <https://github.com/Adrian-Mag/AxiSEM3D_Kernels> (Adrian Mag).
- **Banana-doughnut kernels** — Dahlen, Hung & Nolet (2000); Marquering, Dahlen & Nolet (1999).

Prepared/updated by Jonathan Wolf.
