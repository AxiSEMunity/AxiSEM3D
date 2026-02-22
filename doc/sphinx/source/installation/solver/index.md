# The Solver

The AxiSEM3D solver is built on top of several widely used numerical libraries.
This section describes how to install dependencies and build the solver using CMake.


<span style="color: red;"><b>*** TO DO: This section is a mashup of different documents. Instructions need to be verified and updated to the latest install.</b> </span>
---

## Dependencies

AxiSEM3D depends on the following external libraries.

> **Note:** Only **Eigen** and **Boost** have strict version requirements.  
> Other dependencies (FFTW, METIS, NetCDF, HDF5) are generally not sensitive to minor version differences.

### Required Dependencies

| Name | Role | Version Requirement | Notes |
|------|------|---------------------|------|
| [Eigen](http://eigen.tuxfamily.org/) | Linear algebra | **>= 3.4.0 (strict)** | Current stable release is 3.4 (Mar 2022). |
| [Boost](https://www.boost.org/) | C++ helper utilities | **>= 1.85.0 (strict)** | Only header-only modules are used. |
| [FFTW](http://www.fftw.org/) | Fast Fourier transform | >= 3.3.8 | Version is not sensitive in practice. |
| [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) | Mesh partitioning | >= 5.1.0 | Version is not sensitive; 32/64-bit builds acceptable. |
| [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) | Multi-dimensional I/O | >= 4.4.1 | Parallel NetCDF supported but optional. |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) | Hierarchical data support | Optional | Only required if NetCDF is built with HDF5. |


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
