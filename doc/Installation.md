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

### 1. Dependencies

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
# download Eigen3
git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
# download Boost
wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -x
```

The above lines create a directory `axisem3d_dependencies` that contains `eigen3_develop` and `boost_1_73_0`. 

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the latest stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even one has `Eigen` installed before.

#### 1.2. FFTW, Metis and NetCDF
On a laptop or workstation, they can be easily installed using `conda` (either `miniconda` or `anaconda`):

```bash
conda install -c conda-forge fftw
conda install -c anaconda metis
conda install -c anaconda netcdf4
```

On an HPC cluster, it is most likely that these packages have been installed and optimized due to their popularity. The software packages are usually managed by `module`.

First, check the available packages and versions:
```bash
# check all available packages
module avail
# check available fftw versions
module avail fftw
```

Next, load the required package:
```bash
# load FFTW 3.3.4
module load fftw/3.3.4
```

If a p






Before the installation, one can create a directory to store the dependencies:

Also, if one uses `conda` (either `miniconda` or `anaconda`), its root path can be identified by
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
eyJoaXN0b3J5IjpbNTI4Nzg1NjU0LC0xNTQ5MjI1MjgyLC0xMz
kyNzcwMjE1LDE5NTQ0NTc1MjgsNjUxODMzNjMzLC0xMDgzNTM1
MTAyLDc5MDc0NjM1MSw4Njg3OTY3NDcsNzMzMTcwODI5LC05OT
M5MDU2NzcsLTEzNjEzOTc5MzMsLTIxMTY2NDM4NDIsMTIxNDAy
MTIyLC0xOTMyOTI0Mjc2LC02MzM3NzY5NjQsLTEyNzkzNTQ5MT
QsMTIxNjE5NzE0NSwtMTMyNzAyNjI1MCwtMTM4MTk3NDM2OCw0
NjY4NzA2ODJdfQ==
-->