Adrian Kernel Companion
=======================
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
  │       └── adrian_kernel_companion/   ← you are here
  └── AxiSEM3D_Kernels/      ← axikernels (cloned from Adrian-Mag)

The external package used for the validated run came from:

  https://github.com/Adrian-Mag/AxiSEM3D_Kernels.git

Install `axikernels` before running.  A ready-made conda environment is
provided in `AxiSEM3D_Kernels/`.

From this example directory the axikernels repo is then at:

  ../../../AxiSEM3D_Kernels

Installation (conda — recommended):

  conda env create -f ../../../AxiSEM3D_Kernels/environment.yml
  conda activate axikernels_env
  pip install -e ../../../AxiSEM3D_Kernels

Or with pip only (into any Python ≥ 3.9 environment):

  pip install -e /path/to/AxiSEM3D_Kernels

See the axikernels README for full environment requirements.


Quick start
-----------
1.  Install axikernels (see above).

2.  Copy the compiled axisem3d binary into this folder (or put it on PATH):

        cp /path/to/build/axisem3d .

3.  Run the forward simulation:

        bash run.sh forward

    Output will be written to simu_forward/output/elements/.

    To use more (or fewer) MPI ranks:

        NRANKS=8 bash run.sh forward

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

        kernel_output/vp_kernel.h5    – kernel data (reloadable via axikernels)
        kernel_output/vp_kernel.png   – quick-look figure


Validation status
-----------------
This workflow was re-run successfully against the current upstream AxiSEM3D
build in this workspace in March 2026 using:

  - `bash run.sh forward`
  - `python compute_kernels.py`

Validated revisions in that run:
  - AxiSEM3D example checkout: `30d7ea4ff225db04979ddcb41ada9e6cf16b426b`
  - axikernels checkout: `890084ebe547df8714213b25885810a8be417e9b`

Validated default parameters in that run:
  - forward MPI ranks: `NRANKS=8`
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

The Python side was validated from the `axikernels_env` conda environment
created from `AxiSEM3D_Kernels/environment.yml`.

The solver binary used for validation was the locally built upstream binary
copied from `AxiSEM3D-upstream/build/axisem3d`.  Build that binary using the
normal AxiSEM3D build instructions before running this example.


Figure reproduction
-------------------
The minimum reproducible figure shipped by this example is the quick-look
kernel slice `kernel_output/vp_kernel.png`, generated directly by
`compute_kernels.py`.

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
  - The primary interface is script-based (run.sh / compute_kernels.py).
  - An optional Jupyter notebook demonstrating the same workflow interactively
    lives in AxiSEM3D_Kernels/examples/example_1D_kernel.ipynb.
  - Computational cost for the default setup (50 s period mesh, 8 MPI ranks):
      Forward:  ~12 GB RAM, ~5.8 GB disk
      Backward: ~10 GB RAM, ~4.2 GB disk
      Combined wavefield storage: plan for at least ~10 GB free disk
      Runtime: wall-clock / CPU-hour estimates have not yet been benchmarked
               for this packaged companion example and will vary with MPI
               ranks, CPU generation, and filesystem speed.
  - `bash run.sh forward` refuses to reuse an existing `simu_forward/`
    directory unless you set `CLEAN_FORWARD=1`.


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

