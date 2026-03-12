[![License MIT](https://img.shields.io/badge/License-MIT-red)](https://github.com/AxiSEMunity/AxiSEM3D/blob/main/LICENSE)
[![CI Status](https://github.com/tjhei/AxiSEM3D/actions/workflows/linux.yml/badge.svg)](https://github.com/tjhei/AxiSEM3D/actions/workflows/linux.yml)
[![Docs](https://app.readthedocs.org/projects/axisem3d/badge/)](https://axisem3d.readthedocs.io/)

# AxiSEM3D

AxiSEM3D is a spectral-element based solver for large-scale 3D seismic wave propagation simulations.

See the [documentation](https://axisem3d.readthedocs.io/en/latest/) for details.


## Quick start

Create a new conda environment:

```bash
conda env create -f environment.yml -n axisem3d
```

Activate the environment:

```bash
conda activate axisem3d
```

### Mesher

```bash
pip install "numpy==1.26.4" "scipy==1.13.1" && pip install --no-deps "https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip"
```

### Solver

Configure the build:

```bash
rm -rf build && cmake -B build \
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
