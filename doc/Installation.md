[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation
The installation of `AxiSEM3D` includes three parts: the mesher, the solver and several tools (mostly Python libraries) for pre- and post-processing.

System requirements:
* Unix-like system (`AxiSEM3D` is untested on Windows)
* C++ compiler supporting C++17 (check [C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support))
* Basic development tools: `python`, `pip`, `conda`, `cmake`, `wget`
* MPI (a serial build can be made but is mostly useless) 



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
[`Boost`](https://www.boost.org/) | C++ helpers | 1.73.0 | `AxiSEM3D` only uses some of its header-only modules.
[`FFTW`](http://www.fftw.org/) | fast Fourier transform | 3.3.4 | Both single- and double-precision builds are required.
[`Metis`](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) | mesh partitioning | 5.1.0 | Both 32- and 64-bit builds are acceptable.
[`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/docs/index.html) | efficient multi-dimensional I/O | 4.4.1 | Parallel build is supported but not mandatory.

#### 1.1. Eigen and Boost

`Eigen` and `Boost` are header-only libraries. All one needs to do is to download the source code:
```bash
# create a top working directory
mkdir -p AxiSEM3D_2020
# create dependency directory
mkdir -p AxiSEM3D_2020/dependencies && cd $_
# download Eigen 3.3.9
git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
# download Boost 1.73
wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -jx
```
The above lines will create a directory `AxiSEM3D_2020/dependencies` that contains `eigen3_develop` and `boost_1_73_0`.

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

On an HPC cluster, it is most likely that these packages have been installed and optimized owing to their popularity. On a cluster, the software packages are usually managed by `module`.  To browse available packages and versions:
```bash
# all available packages
module avail
# all available versions of FFTW
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

If a package is missing, one may turn to the admin or install it from scratch following the official instructions ([`FFTW`](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html), [`Metis`](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) and [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)).


### 2. Building AxiSEM3D
#### 2.1. Download the code
```bash
# go to the top working directory
cd $HOME/AxiSEM3D_2020
# download the code
git clone https://github.com/kuangdai/AxiSEM-3D.git AxiSEM3D
```
#### 2.2.  Configure by `cmake`
Before doing `cmake`, one must edit the `_ROOT` variables in `AxiSEM3D/SOLVER/CMakeLists.txt` to point to the correct dependencies, for example, on my own machine (the actual paths are *user-dependent*):

```python
# Eigen and Boost installed by downloading the source code
set(EIGEN3_ROOT $ENV{HOME}/AxiSEM3D_2020/dependencies/eigen3_develop)
set(BOOST_ROOT  $ENV{HOME}/AxiSEM3D_2020/dependencies/boost_1_73_0)
# FFTW, Metis and NetCDF installed by conda
set(FFTW_ROOT   $ENV{HOME}/anaconda3)
set(METIS_ROOT  $ENV{HOME}/anaconda3)
set(NETCDF_ROOT $ENV{HOME}/anaconda3)
```

Alternatively, one can set the corresponding environment variables, leaving `AxiSEM3D/SOLVER/CMakeLists.txt` unchanged:
```bash
# Eigen and Boost installed by downloading the source code
export EIGEN3_ROOT=$HOME/AxiSEM3D_2020/dependencies/eigen3_develop
export BOOST_ROOT=$HOME/AxiSEM3D_2020/dependencies/boost_1_73_0
# FFTW, Metis and NetCDF installed by conda
export FFTW_ROOT=$HOME/anaconda3
export METIS_ROOT=$HOME/anaconda3
export NETCDF_ROOT=$HOME/anaconda3
```
To avoid setting these environment variables every time for a new conversation, one can copy them to `.bash_profile` or `.bashrc`. 

To find the `_ROOT` of a package on an HPC cluster, one can use `module show`, for example:
```bash
# examples for Archer (archer.ac.uk)
# show FFTW
module show fftw/3.3.4.11
# show Metis
module show metis/5.1.0_build2
# show NetCDF
module show cray-netcdf/4.6.1.3
```

<strong>Note</strong>: the `_ROOT` variables sent to `cmake` is neither the library path ended with `/lib` nor the include path ended with `/include`; it is the one containing both `/lib` and `/include`. 

After setting the `_ROOT` variables, one can do `cmake`, sending the C, C++ and Fortran compilers via -Dcc, -Dcxx and -Dftn, respectively: 
```bash
# create a build directory
# It is good practice to put nothing under this 'build'
# directory becuase it must be emptied before redoing 
# cmake when CMakeLists.txt is changed
mkdir -p build && cd $_
# cmake (the build type is Release by default)
# Make sure that the C++ compiler supports C++17.
cmake -Dcc=mpicc -Dcxx=mpicxx -Dftn=mpif90 ../AxiSEM3D/SOLVER
```
Upon a successful `cmake`, a summary will be displayed at the end. Check this summary and make sure that `cmake` has found the correct version of the dependencies. 

#### 2.3.  Compile and link by `make`
To compile and link AxiSEM3D:
```bash
# j8 means using 8 threads to accelerate compilation
make -j8
```

Finally, one can verify the executable:
```bash
# the number of processors can be arbitrary
mpirun -np 4 ./axisem3d
```
That an error message appears saying that "missing input directory"  means AxiSEM3D has been built successfully

If linking fails with missing `H5_`



## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIwODE2NzgzNTMsLTExOTE3MDk3NzIsLT
I5MzgyODE3LC0xNDE4MjAyNzI0LDYwMDYyNDI1MCwxNjE3ODY4
MjI4LC03NjI1MDA2MzksNjEzMzc4ODA1LC0xOTc0MTE0NTcxLC
0xOTExNDQzNzMxLC0yMDQyMjc1MzY1LDE4OTU2MTA3MzksMTkz
NzMyMDk1NywtNDkzNjQ1NTMwLDEzODgxODY0MDIsLTUyMjkxOD
g2MCwtNTQyMTAxMTgzLC0xNjExODM5MDAyLC0xMzE0MjAxNDM5
LC00ODQzOTY3MTRdfQ==
-->