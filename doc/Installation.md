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
This will display all the arguments you can pass to the mesher. 


## Solver

The `AxiSEM3D` solver is developed on top of several modern numerical packages including:

Name|Role|Minimum version|Note
--- | --- | ---|---
`Eigen` | linear algebra | 3.3.9 | The current stable release 3.3.7 (up to July 2020) is insufficient.
`Boost` | C++ helpers | 1.7.1 | `AxiSEM3D` only uses some of its header-only modules.
`FFTW` | fast Fourier transform | 3.3.8 | Both single- and double-precision builds are required.
`Metis` | mesh partitioning | 5.1.0 | Both 32- and 64-bit builds are accepted.
`NetCDF` | intensive structured IO | 4.7.1 | Parallel build is supported but not mandatory.

Before the installation, we can create a directory to store the dependencies:
```bash
$ mkdir -p axisem3d_dependencies && cd $_
```


### 1. Eigen
`Eigen` is a header-only library, so we just need to download the source code:
```bash
$ git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
```
This will create `eigen3_develop` under the current directory. Next, to let `AxiSEM3D` use this version, do
```bash
$ export EIGEN3_ROOT=$PWD/eigen3_develop
``` 

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the current stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even you have installed `Eigen` before. 


### 2. Boost
`AxiSEM3D` only uses some of the header-only modules of `Boost`, so we just need to download the source code:

```bash
$ wget -c https://dl.bintray.com/boostorg/release/1.73.0/source/boost_1_73_0.tar.bz2 -O - | tar -x
```
This will create `eigen3_develop` under the current directory. Next, to let `AxiSEM3D` use this version, do
```bash
$ export EIGEN3_ROOT=$PWD/eigen3_develop
``` 



## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDc2NjQzNjExLC0yMTE2NjQzODQyLDEyMT
QwMjEyMiwtMTkzMjkyNDI3NiwtNjMzNzc2OTY0LC0xMjc5MzU0
OTE0LDEyMTYxOTcxNDUsLTEzMjcwMjYyNTAsLTEzODE5NzQzNj
gsNDY2ODcwNjgyLC0xNjQ3MDc4OTA5LC0xMzgzNzcwMjA2LC0x
NzQ5MDU4NTA1LDEzNzE4ODg1OCwtMzMyNzk0ODY3LC0xNzM3NT
g1MTk1LC01Mjg5MzU5NjEsMTEwNzA2ODY2MCwtMjEwMDQ3MTY0
NywtMjE2MzIxMjM4XX0=
-->