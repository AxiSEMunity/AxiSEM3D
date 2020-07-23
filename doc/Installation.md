[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation
The installation of `AxiSEM3D` includes three parts: the mesher, the solver and several tools (mostly Python libraries) for pre- and post-processing.

System requirements:
* Unix-like system (`AxiSEM3D` is untested on Windows)
* C++ compiler supporting C++17 (check [C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support))
* Basic development tools: `python`, `pip`, `conda` ([`Anaconda`](https://docs.anaconda.com/anaconda/install/) or [`Miniconda`](https://docs.conda.io/en/latest/miniconda.html)), `cmake` (≥ 3.15.3), `wget`
* MPI (a serial build can be made but is not very useful) 



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

The `AxiSEM3D` solver is developed on top of the following modern numerical libraries:

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

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the latest stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are *essential* even one had `Eigen` installed before.


#### 1.2. FFTW, Metis and NetCDF
On a laptop or workstation, these packages can be easily installed using `conda`:

```bash
# install FFTW
conda install -c conda-forge fftw
# install Metis
conda install -c anaconda metis
# install NetCDF
conda install -c anaconda netcdf4
```

On an HPC cluster, it is most likely that these packages have been installed with an optimized configuration owing to their popularity. On a cluster, the software packages are usually managed by `module`.  To browse available packages and versions:
```bash
# all available packages
module avail
# all available versions of FFTW
module avail fftw
```

To load a required package, for example, to load `FFTW` (this is *machine-dependent*):
```bash
module load fftw
```

If a package is missing, one may 
1. turn to admin for help;
2. install it by `conda` (many clusters allow users to install their own `Miniconda` or `Anaconda`);
3. install it from scratch (see instructions for [`FFTW`](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html), [`Metis`](http://glaros.dtc.umn.edu/gkhome/metis/metis/download) and [`NetCDF`](https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html)).  




### 2. Building AxiSEM3D
#### 2.1. Download the code
```bash
# go to the top working directory
cd $HOME/AxiSEM3D_2020
# download the code
git clone https://github.com/kuangdai/AxiSEM-3D.git AxiSEM3D
```

#### 2.2.  Configure by `cmake`
First, create a `build` directory: 
```bash
# under the top working directory
mkdir -p build && cd $_
```
Because `build` must be emptied before redoing `cmake`, it is good practice to put nothing under `build` except those automatically created by `cmake` and `make`. 

Next, do `cmake` like this:
```bash
cmake -Dcc=mpicc -Dcxx=mpicxx -Dftn=mpif90 \
-Deigen=$HOME/AxiSEM3D_2020/dependencies/eigen3_develop \
-Dboost=$HOME/AxiSEM3D_2020/dependencies/boost_1_73_0 \
-Dfftw=$HOME/anaconda3 \
-Dmetis=$HOME/anaconda3 \
-Dnetcdf=$HOME/anaconda3 \
../AxiSEM3D/SOLVER/
```

It can take the following parameters:
Parameter|Role|Default|Note
--- | --- | ---|---
`Dcc`, `Dcxx`, `Dftn`| C, C++, Fortran compilers | gcc, g++, gfortran | The C++ compiler must support C++17
`Deigen`, `Dboost`, `Dfftw`, `Dmetis`, `Dnetcdf`| paths of the dependencies | empty string | Such a path should contain both `\lib` and `\include`. To find the path of a package managed by `module`, use `module show` (e.g., `module show fftw`). 
`Dhdf5` | path of `HDF5` | empty | If `NetCDF` was built as a static library, linking will fail with missing `_H5` symbols. In that case, one has to set `Dhdf5` pointing to the HDF5 library used to build `NetCDF`.





Before running `cmake`, one must edit the `_ROOT` variables in `AxiSEM3D/SOLVER/CMakeLists.txt` to point to the correct dependencies, for example (the actual paths are *user-dependent*):

```python
# Eigen and Boost installed by downloading the source code
set(EIGEN3_ROOT $ENV{HOME}/AxiSEM3D_2020/dependencies/eigen3_develop)
set(BOOST_ROOT  $ENV{HOME}/AxiSEM3D_2020/dependencies/boost_1_73_0)
# FFTW, Metis and NetCDF installed by conda
set(FFTW_ROOT   $ENV{HOME}/anaconda3)
set(METIS_ROOT  $ENV{HOME}/anaconda3)
set(NETCDF_ROOT $ENV{HOME}/anaconda3)
```

Alternatively, one can set the corresponding environment variables, leaving `CMakeLists.txt` unchanged:
```bash
# Eigen and Boost installed by downloading the source code
export EIGEN3_ROOT=$HOME/AxiSEM3D_2020/dependencies/eigen3_develop
export BOOST_ROOT=$HOME/AxiSEM3D_2020/dependencies/boost_1_73_0
# FFTW, Metis and NetCDF installed by conda
export FFTW_ROOT=$HOME/anaconda3
export METIS_ROOT=$HOME/anaconda3
export NETCDF_ROOT=$HOME/anaconda3
```
To avoid setting these environment variables every time in a new conversation, one can copy them to `.bash_profile` or `.bashrc`. 

Use `module show` to find the `_ROOT` of a package (the one containing both `/lib` and `/include`) on a cluster, for example:
```bash
# show FFTW (this is machine-dependent!)
module show fftw
```

After setting the `_ROOT` variables, one can do `cmake`, sending the C, C++ and Fortran compilers via -Dcc, -Dcxx and -Dftn, respectively: 
```bash
# create a build directory
mkdir -p build && cd $_
# cmake (the build type is Release by default)
# make sure that the C++ compiler supports C++17
cmake -Dcc=mpicc -Dcxx=mpicxx -Dftn=mpif90 ../AxiSEM3D/SOLVER
```
Upon a successful `cmake`, a summary will be displayed at the end. Check this summary and make sure that `cmake` has found the correct version of the dependencies. 

<strong>NOTE</strong>: If `NetCDF` was built as a static library, linking (in 2.3) will fail with missing `_H5` symbols. In that case, one has to set `LINK_TO_HDF5` as `true` and provide `HDF5_ROOT`. Also, if `NetCDF` was built statically with remote client support, `-lcurl` must be added to `ADDITIONAL_LIBS`.
 
<strong>NOTE</strong>: Whenever `CMakeLists.txt` has been changed, the build directory must be emptied before redoing `cmake`.

#### 2.3.  Compile and link by `make`
To compile and link AxiSEM3D:
```bash
# -j8 means using 8 threads to accelerate compilation
make -j8
```





Finally, one can verify the executable:
```bash
# the number of processors can be arbitrary
mpirun -np 4 ./axisem3d
```
`AxiSEM3D` has been built successfully if an error message appears saying "Missing input directory".




## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbODAzMjE2MzUwLDM5MzE0NjgyNiwxMTU3OT
AzMzg1LDE1MzY0MzIzNTcsLTE5MjM0NDk2NCwxMjAyMDY4NjIs
Mzg5NDU3MTQ0LDE5NjYwMTQ5OTAsNDMyNzcyMjM4LC0xMjY4Nj
U0NTMsLTc0NTQ0MjUyMiwzMzc2NjIxODUsLTIxODg1MTUyOCwt
MTg3ODk2NzcwMywxMzEwMzc4MzY4LDE5MTI0NTQ5NiwyMDQxND
E4OTkyLDEwODA4NjY3OSwtMTE5MTcwOTc3MiwtMjkzODI4MTdd
fQ==
-->