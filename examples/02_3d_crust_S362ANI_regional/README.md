# Example 02: Regional S362ANI with 1D and 3D crust

This example contains two regional simulations extracted from the same source and mantle setup:

- `simu_with_1d_crust`: mantle heterogeneity from S362ANI with a 1D crust.
- `simu_with_3d_crust`: the same regional setup with a 3D crust derived from Crust1.0.

## Status

Both variants were made runnable with repo-local launchers and validated on this system. The 1D and 3D crust jobs were each submitted on `windfall` and both reached the Newmark time loop.

## Directory layout

- `input_share/`: common mesh and model files for both runs.
- `input_with_1d_crust/`: crust-specific overrides for the 1D case.
- `input_with_3d_crust/`: crust-specific overrides for the 3D case.
- `run_with_1d_crust.sh`, `run_with_3d_crust.sh`: portable launchers.
- `run_with_1d_crust_windfall.sbatch`, `run_with_3d_crust_windfall.sbatch`: Slurm wrappers.
- `gen_crust1.ipynb`: notebook used to build the 3D crust product.
- `post_processing.ipynb`: notebook for comparing the results.

## Prerequisites

1. Build AxiSEM3D at the repository root.
2. Make sure `mpirun` is available.
3. The run scripts assume the input trees already present in this directory.
4. If you are not on a module-based system, remove or adapt the `module load` section in the shell launchers.

## Run locally

From this directory:

```bash
./run_with_1d_crust.sh
./run_with_3d_crust.sh
```

Useful overrides:

```bash
NP=8 ./run_with_1d_crust.sh
NP=8 ./run_with_3d_crust.sh
OUTPUT_DIR=$PWD/simu_with_3d_crust/output_custom ./run_with_3d_crust.sh
```

Each launcher stages a fresh `input/` tree under its simulation directory before starting the solver.

## Run on Slurm

From this directory:

```bash
sbatch run_with_1d_crust_windfall.sbatch
sbatch run_with_3d_crust_windfall.sbatch
```

Each batch wrapper writes to `output_windfall_<jobid>/` inside the matching simulation directory.

## Rebuild the prepared inputs

The regional mesh is already included, but the original mesh command was:

```bash
python -m salvus_mesher_lite.interface AxiSEM \
	--basic.model prem_ani \
	--basic.period 50 \
	--spherical.min_radius 3480 \
	--chunk2D.max_colatitude 60 \
	--output_file regional_mesh__prem_ani__50s.e
```

The S362ANI mantle model was downloaded from IRIS EMC. The 3D crust files in `input_with_3d_crust/` were prepared with `gen_crust1.ipynb`.

## Expected runtime

The original note estimated about 1 hour for the 1D crust run and about 4 hours for the 3D crust run on 4 MPI ranks.

## Post-processing

Use `post_processing.ipynb` to compare the regional seismograms and inspect the impact of the 3D crust.