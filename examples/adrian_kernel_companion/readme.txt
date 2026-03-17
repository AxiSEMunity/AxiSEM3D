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
Install `axikernels` before running.  A ready-made conda environment is
provided in AxiSEM3D_Kernels/:

    conda env create -f AxiSEM3D_Kernels/environment.yml  # creates axikernels_env
    conda activate axikernels_env

Or install into an existing environment:

    pip install -e /path/to/AxiSEM3D_Kernels   # development install

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
  - Computational cost for the default setup (50 s period mesh, 4 MPI ranks):
      Forward:  ~12 GB RAM, ~5.8 GB disk
      Backward: ~10 GB RAM, ~4.2 GB disk
  - `bash run.sh forward` refuses to reuse an existing `simu_forward/`
    directory unless you set `CLEAN_FORWARD=1`.

