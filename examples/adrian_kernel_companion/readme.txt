Adrian Kernel Companion
=======================
This folder scaffolds a companion example for computing 1-D sensitivity
kernels from AxiSEM3D forward and adjoint simulations using the external
Python package `axikernels`.

It is a companion example: AxiSEM3D runs the wave simulations as usual, and
`axikernels` handles the post-processing.  The kernel functionality is NOT
built into AxiSEM3D itself.


Target scope once this example is fully populated
-----------------------------------------------
  - 1-D PREM isotropic elastic model (~40 s period)
  - Single P-wave arrival, Gaussian source-time function
  - Station grid at 10-degree spacing

NOT validated (do not extrapolate this example to):
  - 3-D heterogeneous models
  - Period ranges outside the validated set
  - Kernel types beyond those tested in axikernels


External dependency
-------------------
Before running this example, install axikernels into your Python environment:

  pip install -e /path/to/AxiSEM3D_Kernels   # local development install

See the axikernels README for environment requirements.


Quick start
-----------
This Phase 1 scaffold is not runnable yet.  Phase 2 will populate the real
input files and implement the script-based kernel driver.

Planned user workflow:
1. Copy the compiled axisem3d binary into this folder (or set PATH).
2. Set up and run the forward simulation:
     bash run.sh forward
3. Set up and run the adjoint simulation:
     bash run.sh backward
4. Compute and plot the kernels:
     python compute_kernels.py

Steps 2-3 each create a self-contained simulation directory (simu_forward/,
simu_backward/) and run axisem3d inside it.  Step 4 calls axikernels to
cross-correlate the forward and adjoint wavefields and output kernel slices.


Input files
-----------
  input_forward/   – AxiSEM3D inparam files for the forward run
  input_backward/  – AxiSEM3D inparam files for the adjoint run

These directories are populated in Phase 2 of the setup plan.  For now they
serve as placeholders documenting the expected layout only.


Notes
-----
  - The primary interface is script-based (run.sh / compute_kernels.py).
  - Any optional Jupyter notebook added in Phase 2 will be supplementary only;
    the full workflow should be reproducible without it.
  - This example follows the file-layout conventions of the other examples
    in this directory (run*.sh + readme.txt + input*/).
