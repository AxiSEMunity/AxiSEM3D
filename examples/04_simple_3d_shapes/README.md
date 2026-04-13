# Example 04: Simple 3D shapes

This example collects three small workflows for superimposing 3D heterogeneity on top of a background domain. The notebook `Example_4_README.ipynb` still provides the tutorial narrative, but the shell launchers added here are the fastest way to run the cases on a normal Linux system or on Slurm.

## Status

The following subcases were validated on this system:

- `example_input_cartesian`: launched on `windfall` and reached the Newmark time loop.
- `example_single_plume`: launched on `windfall` and reached the Newmark time loop.
- `example_release_paper`: launcher is working and the current cluster run is still active in a long preloop stage. That case is intentionally being left alone.

## Shared prerequisites

1. Build AxiSEM3D at the repository root.
2. Create and activate a Python environment under `examples/04_simple_3d_shapes/.venv`.
3. Install the Python packages used by the preprocessing scripts:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install 'numpy<2' scipy netCDF4 salvus-mesher-lite
```

`salvus-mesher-lite` currently exposes `python -m salvus_mesher_lite.interface`, not `salvus_mesh_lite`, and it is safer to keep NumPy below 2.0 because the mesher still references `np.infty`.

The run scripts also load site modules if Lmod is present. On systems without those exact module names, remove or adapt that block.

## Subcases

### `example_input_cartesian`

What it does:

- builds a Cartesian background mesh,
- generates a spherical 3D blob with `gen_blob.py`,
- launches the solver with `run.sh` or `run_windfall.sbatch`.

Run locally:

```bash
cd example_input_cartesian
./run.sh
```

Run on Slurm:

```bash
cd example_input_cartesian
sbatch run_windfall.sbatch
```

Compatibility note:

Some `salvus-mesher-lite` builds generate a Cartesian Exodus mesh without the `discontinuities` variable expected by AxiSEM3D. The local `run.sh` checks for that variable and, if needed, falls back to the known-good mesh in `../example_release_paper/input/paper_example.e`.

### `example_single_plume`

What it does:

- generates a global background mesh with `gen_mesh.sh`,
- creates a plume perturbation with `gen_plume.py`,
- optionally prepares a visualization file with `gen_movie.py` after the run.

Run locally:

```bash
cd example_single_plume
./run.sh
```

Run on Slurm:

```bash
cd example_single_plume
sbatch run_windfall.sbatch
```

### `example_release_paper`

What it does:

- builds the paper mesh in `input/`,
- generates a large 3D NetCDF model with `generate_3D_model.py`,
- launches the solver from `run.sh` or `run_windfall.sbatch`,
- reproduces the figure with `reproduce_paper_figure.py` after the run completes.

Run locally:

```bash
cd example_release_paper
./run.sh
```

Run on Slurm:

```bash
cd example_release_paper
sbatch run_windfall.sbatch
```

This is by far the heaviest case in the folder. Expect a long preprocessing and preloop stage, high memory use, and substantially longer queue and runtime requirements than the other two subcases.

## Outputs

Each subcase writes its own `output/` or `output_windfall_<jobid>/` directory inside the subcase folder.

## Tutorial notebook

Use `Example_4_README.ipynb` if you want the step-by-step scientific explanation rather than just the runnable commands.