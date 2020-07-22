[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation
The installation of `AxiSEM3D` includes three parts: the mesher, the solver and several tools (mostly Python libraries) for pre- and post-processing. 


## Mesher
[`SalvusMeshLite`](https://gitlab.com/Salvus/SalvusMeshLite) is the mesher for `AxiSEM3D`. Its installation is trivial with `pip`: 
```bash
$ pip install https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip
```
The mesher runs on a single processor, so there is no need to install it on an HPC cluster. It is efficient enough to generate very large-scale meshes on a laptop.

Verify the installation by
```bash
$ python -m salvus_mesh_lite.interface AxiSEM -h
```
This will display all the arguments one can pass to the mesher. 


## Solver

The `AxiSEM3D` solver is developed on top of several modern numerical packages including:

Name|Role|Minimum version|Note
--- | --- | ---|---
`Eigen` | linear algebra | 3.3.9 | The current stable release 3.3.7 (up to July 2020) is insufficient.
`Boost` | C++ helpers | 1.7.1 | `AxiSEM3D` only uses some of its header-only modules.
`FFTW` | fast Fourier transform | 3.3.8 | Both single- and double-precision builds are required.
`Metis` | mesh partitioning | 5.1.0 | Both 32- and 64-bit builds are accepted.
`NetCDF` | intensive structured I/O | 4.7.1 | Parallel build is supported but not mandatory.

Before the installation, one can create a directory to store the dependencies:
```bash
$ mkdir -p axisem3d_dependencies && cd $_
$ export AXISEM3D_DEPENDS_PATH=$PWD
```
Also, if one uses `conda`, the root path can be identified by
```bash
$ export CONDA_PATH=$(dirname $(dirname $(which conda)))
```

### 1. Eigen
`Eigen` is a header-only library, so one just needs to download the source code:
```bash
$ git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
$ export EIGEN3_ROOT=$AXISEM3D_DEPENDS_DIR/eigen3_develop
```
The first line creates `eigen3_develop` under the current directory. The second line enables `AxiSEM3D` to find this version by setting the environment variable `EIGEN3_ROOT`.

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the current stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even one already has `Eigen` installed before.


### 2. Boost
`AxiSEM3D` only uses some of the header-only modules of `Boost`. Similar to `Eigen`, one only needs to download the source code and set the root variable:

```bash
$ wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -x
$ export BOOST_ROOT=$AXISEM3D_DEPENDS_DIR/boost_1_73_0
```

Alternatively, one can use `conda` to install `Boost`: 
```bash
$ conda install -c conda-forge boost
$ export BOOST_ROOT=$CONDA_PATH
```


### 3. FFTW
Using `conda`:
```bash
$ conda install -c conda-forge fftw
$ export FFTW_ROOT=$CONDA_PATH
```

### 4. Metis
Using `conda`:
```bash
$ conda install -c anaconda metis
$ export METIS_ROOT=$CONDA_PATH
```

### 5. NetCDF
Using `conda`:
```bash
$ conda install -c anaconda netcdf4
$ export NETCDF_ROOT=$CONDA_PATH
$ export HDF5_ROOT=$CONDA_PATH
```
Using a `NetCDF` build with parallel I/O support can enhance the performance of `AxiSEM3D` and simplify post-processing. However, a parallel build is sometimes difficult to make. 


## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNTI3OTYyMDQyLDc5MDc0NjM1MSw4Njg3OT
Y3NDcsNzMzMTcwODI5LC05OTM5MDU2NzcsLTEzNjEzOTc5MzMs
LTIxMTY2NDM4NDIsMTIxNDAyMTIyLC0xOTMyOTI0Mjc2LC02Mz
M3NzY5NjQsLTEyNzkzNTQ5MTQsMTIxNjE5NzE0NSwtMTMyNzAy
NjI1MCwtMTM4MTk3NDM2OCw0NjY4NzA2ODIsLTE2NDcwNzg5MD
ksLTEzODM3NzAyMDYsLTE3NDkwNTg1MDUsMTM3MTg4ODU4LC0z
MzI3OTQ4NjddfQ==
-->