Adrian Kernel 1D
================
A companion example for computing 1-D P-wave sensitivity kernels (banana-
doughnuts) from AxiSEM3D forward simulations using the external Python package
`axikernels`.

This is a companion example: AxiSEM3D runs the forward wave simulation as
usual, and `axikernels` handles adjoint-source construction, backward
simulation, and kernel evaluation.  Kernel computation is NOT built into
AxiSEM3D itself.


Validated scope
---------------
  - 1-D PREM isotropic elastic model
  - Default example uses the 50 s mesh (`prem_iso_elastic_50s.e`)
    with the 40 s mesh available as an alternative input asset
  - Isotropic moment-tensor source (explosive-type)
  - Forward setup includes the 10-degree station grid used in the original
    axikernels example assets
  - P400P phase cross-correlation time-shift kernel (vp)

NOT validated (do not extrapolate this example to):
  - 3-D heterogeneous models
  - Period ranges outside the tested set
  - Kernel types beyond those listed in the axikernels README


External dependency
-------------------
This companion example assumes that `AxiSEM3D_Kernels` has been cloned as a
separate repository.  Clone both repos side by side so you end up with:

  some_folder/
  ├── AxiSEM3D/              ← this repository (cloned from AxiSEMunity)
  │   └── examples/
  │       └── 07_kernels/
  │           └── kernel_1D/   ← you are here
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

2.  Build AxiSEM3D and copy the compiled axisem3d binary into this folder
    (or put it on PATH):

        cp /path/to/axisem3d .

3.  Run the forward simulation.  axikernels needs one simulation directory
    holding input/, output/, and the axisem3d binary (which it reuses for the
    backward run), so assemble simu_forward/ and run the solver inside it:

        mkdir -p simu_forward
        cp -r input_forward simu_forward/input
        cp /path/to/axisem3d simu_forward/axisem3d
        mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output

    Output will be written to simu_forward/output/elements/.  Adjust -np to the
    number of MPI ranks you want.

4.  Compute the vp sensitivity kernel:

        python compute_kernels.py

    This step:
      - constructs the adjoint source automatically from the forward wavefield
      - creates and runs the backward simulation (no manual invocation needed)
      - evaluates the vp kernel on a 2-D great-circle slice
      - writes kernel_output/vp_kernel.h5 and kernel_output/vp_kernel.png
      - auto-confirms the current interactive axikernels prompts

    Re-running this step will overwrite the generated backward simulation
    directory `backward_simu_forward/` unless you rename or move it first.

    Run with --help to see adjustable parameters:

        python compute_kernels.py --help

5.  Inspect outputs:

        kernel_output/vp_kernel.h5    – vp kernel data (reloadable via axikernels)
        kernel_output/vp_kernel.png   – quick-look vp kernel figure
        kernel_output/moho_kd.h5      – Moho interface kernel data (see below)
        kernel_output/moho_kd.png     – quick-look Moho kernel figure


Validation status
-----------------
This workflow was validated end to end against a recent upstream AxiSEM3D
build using:

  - mpirun -np 8 axisem3d --input simu_forward/input --output simu_forward/output
  - python compute_kernels.py

Validated default parameters in that run:
  - forward MPI ranks: `-np 8`
  - backward MPI ranks: `--cores 8`
  - mesh: `prem_iso_elastic_50s.e`
  - receiver: `--receiver 0 40`
  - phase window: `--window 425 475`
  - time shift: `--tau 2`
  - channel: `--channel UZ`
  - slice resolution: `--resolution 200`

The validated run:
  - completed the forward simulation to 600.223 s
  - generated and ran `backward_simu_forward/` automatically
  - wrote `kernel_output/vp_kernel.h5`
  - wrote `kernel_output/vp_kernel.png`

The Python side was validated from the axikernels environment described in
`AxiSEM3D_Kernels/environment.yml`.

Build the AxiSEM3D solver binary using the normal AxiSEM3D build instructions
before running this example.


Moho topography sensitivity kernel (K_SS)
------------------------------------------
`compute_kernels.py` also computes a Moho interface sensitivity kernel after
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

  moho_kd.png  – quick-look pcolormesh figure of K_SS with source (*) and
                 receiver (^) positions marked; colourscale is symmetric RdBu_r
                 with the range set to 10 % of the absolute maximum.


Figure reproduction
-------------------
The minimum reproducible figures shipped by this example are the quick-look
kernel slice `kernel_output/vp_kernel.png` and the Moho interface kernel map
`kernel_output/moho_kd.png`, both generated directly by `compute_kernels.py`.

This is intended to demonstrate the validated companion workflow end to end.
It is not yet a full packaging of every publication-quality figure previously
produced during the original axikernels work.  The interactive notebook in
`AxiSEM3D_Kernels/examples/example_1D_kernel.ipynb` remains useful for more
exploratory visualization and figure customization.


No manual backward run required
--------------------------------
The backward (adjoint) simulation is created and launched automatically by
`axikernels`.  You do not need to populate `input_backward/` or call any
separate run script.  The directory `simu_forward/` contains everything
axikernels needs (input files + axisem3d binary) to mirror the setup.

If you pass a custom `--forward` path to `compute_kernels.py`, it must point
to the full simulation directory, not only to `output/elements/`.


Input files
-----------
  input_forward/   – AxiSEM3D inparam files for the 1-D PREM forward run
                     (populated from the validated axikernels 1D kernel example)
  input_backward/  – NOT NEEDED; see input_backward/PLACEHOLDER.txt for details


Notes
-----
  - The forward run is a plain mpirun invocation; compute_kernels.py drives the
    adjoint-source construction, backward run, and kernel evaluation.
  - An optional Jupyter notebook demonstrating the same workflow interactively
    lives in AxiSEM3D_Kernels/examples/example_1D_kernel.ipynb.
  - Computational cost for the default setup (50 s period mesh, 8 MPI ranks):
      Forward:  ~12 GB RAM, ~5.8 GB disk
      Backward: ~10 GB RAM, ~4.2 GB disk
      Combined wavefield storage: plan for at least ~10 GB free disk
      Runtime: wall-clock / CPU-hour estimates have not yet been benchmarked
               for this packaged companion example and will vary with MPI
               ranks, CPU generation, and filesystem speed.
  - To re-run from scratch, delete simu_forward/, backward_simu_forward/, and
    kernel_output/ first (all git-ignored).


Integration note
----------------
This example is deliberately packaged as a companion example, not as evidence
that kernel computation is integrated into AxiSEM3D core.

What is demonstrated here:
  - AxiSEM3D can provide the forward and adjoint wavefields required for a
    kernel workflow.
  - An external package (`axikernels`) can turn those wavefields into a
    concrete sensitivity-kernel example.

What remains outside the scope of this example:
  - General in-core kernel support in AxiSEM3D itself
  - Broad coverage of kernel types, objectives, and model classes
  - Full replacement of the exploratory notebook-based workflow

