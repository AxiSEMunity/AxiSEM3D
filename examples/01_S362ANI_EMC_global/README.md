# Example 01: Global 3D S362ANI

This example adds the S362ANI mantle model to a global PREM background and uses the 2011 Virginia earthquake plus GSN and USArray stations. The mesh and the 3D model file are already staged in `input/`.

## Status

This case was validated on this system with the repo-local launcher. It starts cleanly locally and was also submitted on Slurm to the `windfall` partition, where it reached the Newmark time loop.

## Directory layout

- `input/`: ready-to-run mesh, 3D model, source, stations, and YAML files.
- `run.sh`: portable launcher that uses `../../build/axisem3d` instead of a copied binary.
- `run_windfall.sbatch`: Slurm wrapper used on this system.
- `output/`, `output_smoke/`, `output_windfall_304982/`: local and cluster outputs.
- `post_processing.ipynb`: notebook for seismograms and animations.

## Prerequisites

1. Build AxiSEM3D at the repository root.
2. Make sure `mpirun` is available.
3. If you use the provided scripts on another cluster, either keep the module names aligned with your site or remove that module-loading block and preload your environment some other way.

## Run locally

From this directory:

```bash
./run.sh
```

Useful overrides:

```bash
NP=8 ./run.sh
AXISEM3D_BIN=/path/to/axisem3d ./run.sh
OUTPUT_DIR=$PWD/output_custom ./run.sh
```

## Run on Slurm

From this directory:

```bash
sbatch run_windfall.sbatch
```

The batch file defaults to the repo build and writes output to `output_windfall_<jobid>/`.

## Rebuild the mesh or 3D model

The prepared inputs are already present, but the original setup was:

```bash
python -m salvus_mesher_lite.interface AxiSEM \
	--basic.model prem_ani \
	--basic.period 50 \
	--output_file input/global_mesh__prem_ani__50s.e
```

The S362ANI NetCDF file in `input/` came from IRIS EMC.

## Expected runtime

The original note estimated about 35 minutes on 4 MPI ranks. Cluster runtime depends strongly on core count and filesystem speed.

## Post-processing

Use `post_processing.ipynb` to inspect seismograms and compare this 3D run against the 1D reference in example 00.