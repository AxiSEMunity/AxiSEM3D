Adrian Kernel 3D
================
A companion example for computing P-wave sensitivity kernels (banana-
doughnuts) and a Moho topography interface kernel from an AxiSEM3D **3D**
forward simulation using the external Python package `axikernels`.

The 3D simulation uses:
  - S362ANI volumetric seismic heterogeneity (×2 amplitude exaggeration)
  - A Gaussian Moho topography perturbation (±5 km, StructuredGridG3D)

This is a companion example: AxiSEM3D runs the forward wave simulation as
usual, and `axikernels` handles adjoint-source construction, backward
simulation, and kernel evaluation.  Kernel computation is NOT built into
AxiSEM3D itself.


Validated scope
---------------
  - 3-D model: S362ANI heterogeneity + Gaussian Moho topography
  - 50 s mesh (`prem_iso_elastic_50s.e`) with PREM isotropic elastic
    background
  - Isotropic moment-tensor source (explosive-type)
  - Forward setup includes the 10-degree station grid
  - P400P phase cross-correlation time-shift kernel (vp)
  - Moho topography interface kernel (K_SS) with topography contour overlay

External dependency
-------------------
This companion example assumes that `AxiSEM3D_Kernels` has been cloned as a
separate repository.  Clone both repos side by side so you end up with:

  some_folder/
  ├── AxiSEM3D/              ← this repository (cloned from AxiSEMunity)
  │   └── examples/
  │       └── 07_kernels/
  │           └── kernel_3D/   ← you are here
  └── AxiSEM3D_Kernels/      ← axikernels (cloned from Adrian-Mag)

The external package used for the validated run came from:

  https://github.com/Adrian-Mag/AxiSEM3D_Kernels.git

Install `axikernels` before running.  A ready-made conda environment is
provided in `AxiSEM3D_Kernels/`.

From this example directory the axikernels repo is then at:

  ../../../../AxiSEM3D_Kernels

Installation (conda — recommended):

  conda env create -f ../../../../AxiSEM3D_Kernels/environment.yml
  conda activate <env-name>          # use the name defined in that environment.yml
  pip install -e ../../../../AxiSEM3D_Kernels

Or with pip only (into any Python ≥ 3.9 environment):

  pip install -e /path/to/AxiSEM3D_Kernels

See the axikernels README for full environment requirements.


Quick start
-----------
1.  Install axikernels (see above).

2.  Provide the NetCDF input files.  Two model files are required in
    input_forward/:

        input_forward/S362ANI_radius.nc      – S362ANI model on radius grid
        input_forward/moho_topography.nc     – Gaussian Moho topography

    moho_topography.nc ships.  S362ANI_radius.nc (~9.4 MB) does NOT ship (it
    exceeds the size limit); regenerate it (or download it from the gallery)
    before running.  To regenerate either file from scratch:

        python convert_s362ani_to_radius.py     # -> input_forward/S362ANI_radius.nc
        python create_moho_topography.py        # -> input_forward/moho_topography.nc

3.  Build AxiSEM3D and copy the compiled axisem3d binary into this folder
    (or put it on PATH):

        cp /path/to/axisem3d .

4.  Run the forward simulation.  axikernels needs one simulation directory
    holding input/, output/, and the axisem3d binary (which it reuses for the
    backward run), so assemble simu_forward/ and run the solver inside it:

        mkdir -p simu_forward
        cp -r input_forward simu_forward/input
        cp /path/to/axisem3d simu_forward/axisem3d
        mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output

    Output will be written to simu_forward/output/elements/.  Adjust -np to the
    number of MPI ranks you want.

5.  Compute the vp sensitivity kernel and Moho interface kernel:

        python compute_kernels_3D.py

    This step:
      - constructs the adjoint source automatically from the forward wavefield
      - creates and runs the backward simulation (no manual invocation needed)
      - evaluates the vp kernel on a 2-D great-circle slice
      - evaluates the Moho interface kernel (K_SS) with topography contours
      - writes all outputs to kernel_output/

    Re-running will overwrite the generated backward simulation directory
    `backward_simu_forward/` unless you rename or move it first.

    Run with --help to see adjustable parameters:

        python compute_kernels_3D.py --help

6.  Inspect outputs:

        kernel_output/vp_kernel.h5    – vp kernel data (reloadable via axikernels)
        kernel_output/vp_kernel.png   – quick-look vp kernel figure
        kernel_output/moho_kd.h5      – Moho interface kernel data (see below)
        kernel_output/moho_kd.png     – Moho kernel figure with topography overlay


Validation status
-----------------
This workflow was validated end to end against a recent upstream AxiSEM3D
build using:

  - mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output
  - python compute_kernels_3D.py

Validated default parameters in that run:
  - forward MPI ranks: `-np 8`
  - backward MPI ranks: `--cores 8`
  - mesh: `prem_iso_elastic_50s.e`
  - receiver: `--receiver 0 40`
  - phase window: `--window 425 475`
  - time shift: `--tau 2`
  - channel: `--channel UZ`
  - slice resolution: `--resolution 200`

The Python side was validated from the axikernels environment described in
`AxiSEM3D_Kernels/environment.yml`.

Build the AxiSEM3D solver binary using the normal AxiSEM3D build instructions
before running this example.


Moho topography sensitivity kernel (K_SS)
------------------------------------------
`compute_kernels_3D.py` computes a Moho interface sensitivity kernel after
the vp kernel, reusing the same forward and backward wavefields with no
additional simulation runs.

The Moho is modelled as a solid-solid (SS) discontinuity at the PREM radius
r = 6 346 600 m (depth 24.4 km), confirmed from the `prem_iso_elastic.bm`
file shipped with the example.

The script evaluates three complementary quantities on a regular lat/lon grid
covering the source–receiver great-circle path:

  K_SS   – total interface kernel (normal + volumetric)
           This is the net sensitivity to a perturbation of the Moho depth:
           a non-zero K_SS at a surface point means that horizontally
           displacing the Moho boundary at that latitude/longitude changes
           the selected travel-time measurement (P400P cross-correlation
           time shift).  Negative values indicate that a locally deeper Moho
           reduces the measured time shift; positive values the opposite.

  K_dn   – normal-displacement (traction-jump) part of the interface kernel.
           This term comes from the jump in traction across the boundary and
           captures sensitivity to the orientation of the boundary normal.

  K_dv   – velocity-contrast (material-jump) part of the interface kernel.
           This term accounts for the contrast in elastic parameters (Vp, Vs,
           rho) across the Moho and would vanish if the two sides were
           identical.

The Moho kernel figure (moho_kd.png) additionally overlays contour lines of
the Gaussian Moho topography from `input_forward/moho_topography.nc`,
providing visual confirmation that the boundary kernel correlates with the
undulation pattern.

Outputs written to `kernel_output/`:

  moho_kd.h5   – HDF5 file containing the 2-D lat/lon grids and the three
                 kernel component arrays:
                   lat_deg  – latitude meshgrid (degrees)
                   lon_deg  – longitude meshgrid (degrees)
                   K_SS     – total interface kernel
                   K_dn     – normal-displacement part
                   K_dv     – velocity-contrast part
                 Attributes:
                   moho_radius_m – 6346600.0 (PREM Moho radius in metres)

  moho_kd.png  – quick-look pcolormesh figure of K_SS with topography contours,
                 source (*) and receiver (^) marked; colourscale is symmetric
                 RdBu_r with range set to 10 % of the absolute maximum.


Input files
-----------
  input_forward/                  – AxiSEM3D inparam files for the 3-D forward run
    inparam.model.yaml            – References S362ANI_radius.nc and
                                    moho_topography.nc for the 3-D model
    S362ANI_radius.nc             – S362ANI volumetric heterogeneity (×2 exag.)
                                    generated by convert_s362ani_to_radius.py
    moho_topography.nc            – Gaussian Moho topography (±5 km)
                                    generated by create_moho_topography.py
    prem_iso_elastic_50s.e        – 50 s PREM isotropic mesh
    prem_iso_elastic.bm           – PREM background model
    STA_10DEG_GRID.txt            – 10-degree station grid
    inparam.advanced.yaml
    inparam.nr.yaml
    inparam.output.yaml
    inparam.source.yaml


Notes on coordinate system
---------------------------
AxiSEM3D stores element output in **reference** coordinates even for 3-D runs.
The PRT transformation (which maps the physical 3-D geometry to the
axisymmetric reference mesh) modifies the physics internally; coordinates
in the output are always referenced to the undeformed background mesh.
Therefore `MOHO_RADIUS` is the correct evaluation radius in both 1-D and 3-D
cases, and the same axikernels kernel evaluation code works for both.


No manual backward run required
---------------------------------
The backward (adjoint) simulation is created and launched automatically by
`axikernels` when you run `compute_kernels_3D.py`.  You do not need to run any
separate backward simulation script.


Notes
-----
  - Computational cost for the default setup (50 s period mesh, 8 MPI ranks):
      Forward:  ~12–14 GB RAM, ~6–8 GB disk (3-D run is heavier than 1-D)
      Backward: ~10 GB RAM, ~4–5 GB disk
      Combined wavefield storage: plan for at least ~14 GB free disk
      Runtime: will vary with MPI ranks, CPU generation, and filesystem speed;
               3-D runs are significantly more expensive than the 1-D example.
  - To re-run from scratch, delete simu_forward/, backward_simu_forward/, and
    kernel_output/ first (all git-ignored).
  - The preprocessing scripts (`convert_s362ani_to_radius.py`,
    `create_moho_topography.py`) write their outputs directly to
    `input_forward/`.  The committed NC files in this directory are the
    outputs from those scripts; you normally do not need to re-run them.


Integration note
----------------
This example is deliberately packaged as a companion example, not as evidence
that kernel computation is integrated into AxiSEM3D core.

What is demonstrated here:
  - AxiSEM3D can run 3-D forward simulations with volumetric heterogeneity
    and Moho topography.
  - An external package (`axikernels`) can turn those wavefields into
    concrete sensitivity kernels (volumetric vp and Moho interface).

What remains outside the scope of this example:
  - General in-core kernel support in AxiSEM3D itself
  - Broad coverage of kernel types, objectives, and model classes
  - Full replacement of the exploratory notebook-based workflow
