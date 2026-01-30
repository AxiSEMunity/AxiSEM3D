# AxiSEM3D

AxiSEM3D is a spectral-element based solver for large-scale 3D seismic wave propagation simulations.

## Quick start

Create a new conda environment:

```bash
conda env create -f environment.yml -n axisem3d
```

Activate the environment:

```bash
conda activate axisem3d
```

Configure the build:

```bash
rm -rf build && cmake -S ./SOLVER -B build \
  -Dcxx=mpicxx \
  -Dhdf5=$CONDA_PREFIX \
  -Dnetcdf=$CONDA_PREFIX \
  -Deigen=$CONDA_PREFIX \
  -Dboost=$CONDA_PREFIX \
  -Dfftw=$CONDA_PREFIX \
  -Dmetis=$CONDA_PREFIX
```

Compile and link:

```bash
cmake --build build -j4
```

Check the executable:

```bash
./build/axisem3d --help
```

For installation on HPC clusters, please refer to the user guide and `tools/installation_scripts`.
