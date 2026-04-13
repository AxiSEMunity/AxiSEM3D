# Example 00: Global 1D PREM

This example runs a 50 s global simulation using the 1D anisotropic PREM background model for the 2011 Virginia earthquake. It is the simplest end-to-end example in the repository and a good first validation target after building AxiSEM3D.

## Status

The completed reference output is stored in `output_test/`.

## Directory layout

- `input/`: prepared mesh, source, stations, and runtime YAML files.
- `output_test/`: completed reference output from the validated run.
- `post_processing.ipynb`: notebook for plotting seismograms and animations.

## Prerequisites

1. Build AxiSEM3D at the repository root so `build/axisem3d` exists.
2. Make sure `mpirun` works in your environment.
3. If you want to regenerate the mesh instead of using the provided one, install `salvus-mesher-lite` and run the mesher command shown below.

On module-based clusters, load your compiler, MPI, NetCDF, FFTW, and METIS stack before running. On systems without modules, just make sure those libraries are already visible to the built binary.

## Run locally

From this directory:

```bash
mpirun -np 4 ../../build/axisem3d \
	--input input \
	--output output
```

Change `-np 4` to match the MPI size you want.

## Run on Slurm

There is no dedicated batch wrapper in this directory, but the case can be submitted directly:

```bash
sbatch --partition=<partition> --nodes=1 --ntasks=4 --time=01:00:00 \
	--wrap='cd <repo>/examples/00_global_1D && mpirun -np ${SLURM_NTASKS} ../../build/axisem3d --input input --output output_${SLURM_JOB_ID}'
```

## Regenerate the mesh

The mesh in `input/` is already usable. If you want to rebuild it:

```bash
python -m salvus_mesher_lite.interface AxiSEM \
	--basic.model prem_ani \
	--basic.period 50 \
	--output_file input/global_mesh__prem_ani__50s.e
```

## Expected output

The run creates an output directory containing station data, diagnostics, and optional wavefield products requested by the YAML files in `input/`.

## Post-processing

Open `post_processing.ipynb` to inspect the seismograms and map-based visualizations for this 1D reference simulation.