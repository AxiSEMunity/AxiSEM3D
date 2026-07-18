[![License MIT](https://img.shields.io/badge/License-MIT-red)](https://github.com/AxiSEMunity/AxiSEM3D/blob/main/LICENSE)
[![CI Status](https://github.com/tjhei/AxiSEM3D/actions/workflows/linux.yml/badge.svg)](https://github.com/tjhei/AxiSEM3D/actions/workflows/linux.yml)
[![Docs](https://app.readthedocs.org/projects/axisem3d/badge/)](https://axisem3d.readthedocs.io/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20142088.svg)](https://doi.org/10.5281/zenodo.20142088)

# AxiSEM3D

AxiSEM3D is a spectral-element based solver for large-scale 3D seismic wave propagation simulations.

See the [documentation](https://axisem3d.readthedocs.io/en/latest/) for details.

Our website is under development and can be found [here](https://axisemunity.github.io/AxiSEM3D/).            

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
rm -rf build && cmake -B build
```

Compile and link:

```bash
cmake --build build -j4
```

Check the executable:

```bash
./build/axisem3d --help
```

## Integration tests

AxiSEM3D comes with a testsuite of integration tests that run specific configurations specified in ``./tests/`` and compare them to reference data. For this, the tools ``numdiff`` is needed.

On x86 machines you can add it to your conda environment using
```bash
conda install -n axisem3d conda-forge::numdiff
```
but this package currently does not exist for ARM (Apple machines). You can use homebrew instead:
```bash
brew install numdiff
```

To run the integration tests:

```bash
ctest --test-dir build --output-on-failure
```

For installation on HPC clusters, please refer to the user guide and `tools/installation_scripts`.
