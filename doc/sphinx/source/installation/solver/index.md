# The Solver

The AxiSEM3D solver is built on top of several widely used numerical libraries.
This section describes how to install these dependencies and how to configure and build the project using CMake.

## Dependencies

AxiSEM3D depends on the following external libraries that need to be installed in order to configure AxiSEM3D. We will describe several ways how to obtain these dependencies below.

> **Note:** Only **Eigen** and **Boost** have strict version requirements that are checked during
configuration. All other dependencies (FFTW, METIS, NetCDF, HDF5) are generally not sensitive to minor version differences
and the versions indicated below are suggestions that are known to work.

| Name | Role | Version Requirement | Notes |
|------|------|---------------------|------|
| [Eigen](http://eigen.tuxfamily.org/) | Linear algebra | **>= 3.4.0 (strict)** | Suggested version 3.4.1 (September 2025). |
| [Boost](https://www.boost.org/) | C++ helper utilities | **>= 1.83.0 (strict)** | Only header-only modules are used. |
| [FFTW](http://www.fftw.org/) | Fast Fourier transform | 3.3.8 | Version is not sensitive in practice. |
| [METIS](https://github.com/KarypisLab/METIS) | Mesh partitioning | 5.1.0 | Version is not sensitive; 32/64-bit builds acceptable. |
| [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) | Multi-dimensional I/O | 4.4.1 | Parallel NetCDF supported but optional. |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) | Hierarchical data support | Optional | Only required if NetCDF is built with HDF5. |


## CMake configuration

There are several ways to install the above-mentioned dependencies. On laptops and workstations a [local installation](./local.md) using Conda is probably the
easiest way. If necessary, the following CMake options can be specified (using the `-D OPTION=VALUE` syntax):

| Variable | Default | Description |
|----------|---------|-------------|
| ADDITIONAL_CXX_FLAGS | | Optional compilation flags (for example for optimization) |
| ADDITIONAL_LINK_OPTIONS | | Optional additional link options |
| SERIAL_BUILD | OFF | Turn to ON to disable MPI |
| NPOL_SEM | 4 | Polynomial order for spectral elements (from 1 to 8) |
| USE_DOUBLE | OFF | If ON use double precision, otherwise single precision numbers |
| MEMORY_SAVING_MODE | OFF | If ON, try to conserve memory |
| SKIP_MM_SET | OFF | If ON, skip _MM_SET_FLUSH_ZERO_MODE in fenv.cpp |
| BOOST_DIR | | Optional hint to the Boost directory |
| EIGEN3_DIR | | Optional hint to the Eigen directory |
| FFTW3_DIR | | Optional hint to the FFTW directory |
| NETCDF_DIR | | Optional hint to the NetCDF directory |
| USE_PARALLEL_NETCDF | OFF | If set to ON, AxiSEM3D will write netCDF output with parallel routines. |
| WITH_HDF5 | ON | If set to OFF or HDF5 is not found, disable HDF5 support |
| HDF5_DIR | | Optional hint to the HDF5 directory |
| COMPILE_DEBUG_EXECUTABLE | ON | Set to OFF to stop cmake from compilining axisem3d-debug executable alongside the optimized axisem3d executable |
---

```{toctree}
---
maxdepth: 1
---
local.md
hpc.md
docker.md
example-install-cases.md
```
