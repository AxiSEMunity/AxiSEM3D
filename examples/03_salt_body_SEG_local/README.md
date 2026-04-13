# Example 03: SEG Salt Body, Cartesian local model

This example contains two local Cartesian simulations based on the SEG C3 salt model:

- `simu1D`: background-only reference case.
- `simu3D`: 3D salt-body perturbation case.

## Status

Both variants were updated to use the repo-local solver and both were validated on this system. The 1D and 3D jobs were submitted on `windfall` and both reached the Newmark time loop.

## Directory layout

- `input1D/`: complete staged input for the 1D case.
- `input3D/`: complete staged input for the 3D case, including `SEG_C3_data.tar.bz2`.
- `run1D.sh`, `run3D.sh`: local launchers.
- `run1D_windfall.sbatch`, `run3D_windfall.sbatch`: Slurm wrappers.
- `SEG_C3.bm`: background model used to generate the Cartesian mesh.
- `post_processing.ipynb`: notebook for visualizing the result.

## Prerequisites

1. Build AxiSEM3D at the repository root.
2. Make sure `mpirun` is available.
3. If you regenerate the mesh, install `salvus-mesher-lite`.

The 3D launcher deliberately extracts `SEG_C3_data.tar.bz2` with Python's `tarfile` module so it does not depend on `lbzip2` or a specific `bzip2` path.

## Run locally

From this directory:

```bash
./run1D.sh
./run3D.sh
```

Useful overrides:

```bash
NP=8 ./run3D.sh
OUTPUT_DIR=$PWD/simu3D/output_custom ./run3D.sh
```

## Run on Slurm

From this directory:

```bash
sbatch run1D_windfall.sbatch
sbatch run3D_windfall.sbatch
```

## Rebuild the mesh

The prepared inputs are already present, but the original mesh command was:

```bash
python -m salvus_mesher_lite.interface AxiSEMCartesian \
	--basic.model SEG_C3.bm \
	--basic.period .2 \
	--cartesian2Daxisem.x 5. \
	--cartesian2Daxisem.min_z 6367.0 \
	--attenuation.frequencies 0.01 10. \
	--output_filename local_mesh__SEG_salt__5Hz.e
```

## Expected runtime

The original note estimated about 0.3 minutes for the 1D case and about 10 minutes for the 3D case on 4 MPI ranks.

## Post-processing

Use `post_processing.ipynb` to inspect the local wavefield response and compare the 1D and 3D runs.