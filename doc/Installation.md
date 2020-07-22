[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation
The installation of `AxiSEM3D` includes three parts: the mesher, the solver and several tools (mostly Python libraries) for pre- and post-processing.


## Mesher
[`SalvusMeshLite`](https://gitlab.com/Salvus/SalvusMeshLite) is the mesher for `AxiSEM3D`. Its installation is trivial with `pip`: 
```bash
pip install https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip
```
The mesher runs on a single processor, so there is no need to install it on an HPC cluster. It is efficient enough to generate very large-scale meshes on a laptop.

Verify the installation by
```bash
python -m salvus_mesh_lite.interface AxiSEM -h
```
This will display all the arguments one can pass to the mesher. 


## Solver

### 1. Installing dependencies

The `AxiSEM3D` solver is developed on top of several modern numerical packages including:

Name|Role|Minimum version|Note
--- | --- | ---|---
[`Eigen`](http://eigen.tuxfamily.org/index.php?title=Main_Page) | linear algebra | 3.3.9 | The current stable release 3.3.7 (up to July 2020) is insufficient.
[`Boost`](https://www.boost.org/) | C++ helpers | 1.7.1 | `AxiSEM3D` only uses some of its header-only modules.
[`FFTW`](http://www.fftw.org/) | fast Fourier transform | 3.3.4 | Both single- and double-precision builds are required.
[`Metis`](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) | mesh partitioning | 5.1.0 | Both 32- and 64-bit builds are acceptable.
[`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/docs/index.html) | efficient multi-dimensional I/O | 4.4.1 | Parallel build is supported but not mandatory.

#### 1.1. Eigen and Boost

`Eigen` and `Boost` are header-only libraries. All one needs to do is to download the source code:
```bash
# create a directory to store Eigen and Boost
mkdir -p axisem3d_dependencies && cd $_
# download Eigen 3.3.9
git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
# download Boost 1.73
wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -jx
```
The above lines will create a directory `axisem3d_dependencies` that contains `eigen3_develop` and `boost_1_73_0`.

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the latest stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are *essential* even one has had `Eigen` installed before.


#### 1.2. FFTW, Metis and NetCDF
On a laptop or workstation, these packages can be easily installed using `conda` (either [`Anaconda`](https://docs.anaconda.com/anaconda/install/) or [`Miniconda`](https://docs.conda.io/en/latest/miniconda.html)):

```bash
# install FFTW
conda install -c conda-forge fftw
# install Metis
conda install -c anaconda metis
# install NetCDF
conda install -c anaconda netcdf4
```

On an HPC cluster, it is most likely that these packages have been installed and optimized owing to their popularity. On a cluster, the software packages are usually managed by `module`.  To check available packages and versions:
```bash
# check available packages
module avail
# check available versions of FFTW
module avail fftw
```

To load the required packages (note that the following lines are *machine-dependent*):
```bash
# examples for Archer (archer.ac.uk)
# load FFTW
module load fftw/3.3.4.11
# load Metis
module load metis/5.1.0_build2
# load NetCDF
module load cray-netcdf/4.6.1.3
```

If one of them is missing, you may turn to the admin for help or install it from scratch following the official instructions ([`FFTW`](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html), [`Metis`](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) and [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)).


### 2. Build AxiSEM3D



```bash
$ export CONDA_PATH=$(dirname $(dirname $(which conda)))
```

### 1. Eigen
`Eigen` is a header-only library, so one just needs to download the source code:
```bash
$ git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
$ export EIGEN3_ROOT=$AXISEM3D_DEPENDS_PATH/eigen3_develop
```
The first line creates `eigen3_develop` under the current directory. The second line enables `cmake` to find this version by setting the environment variable `EIGEN3_ROOT`.

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the latest stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even one already has `Eigen` installed before.


### 2. Boost
`AxiSEM3D` only uses some of the header-only modules of `Boost`. Similar to `Eigen`, one only needs to download the source code and set the root variable:

```bash
$ wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -x
$ export BOOST_ROOT=$AXISEM3D_DEPENDS_PATH/boost_1_73_0
```

Alternatively, one can use `conda`: 
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
Using a `NetCDF` build with parallel I/O support can enhance the performance of `AxiSEM3D` and simplify post-processing. However, a parallel build is sometimes difficult to make. Instructions are provided [here](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html#build_parallel). 


## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTMzMTQ0MTYyMCwtMTYxMTgzOTAwMiwtMT
MxNDIwMTQzOSwtNDg0Mzk2NzE0LDEyNTU0MjI5NjQsLTYyMTY1
ODgxNCwtMTU0OTIyNTI4MiwtMTM5Mjc3MDIxNSwxOTU0NDU3NT
I4LDY1MTgzMzYzMywtMTA4MzUzNTEwMiw3OTA3NDYzNTEsODY4
Nzk2NzQ3LDczMzE3MDgyOSwtOTkzOTA1Njc3LC0xMzYxMzk3OT
MzLC0yMTE2NjQzODQyLDEyMTQwMjEyMiwtMTkzMjkyNDI3Niwt
NjMzNzc2OTY0XX0=
-->