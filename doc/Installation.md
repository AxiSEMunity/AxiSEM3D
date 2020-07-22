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

which will display all the arguments you can pass to the mesher. 


## Solver

The `AxiSEM3D` solver is developed on top of several modern numerical packages including `Eigen`, `Boost`,  `FFTW`, `Metis` and `NetCDF`. 

### 1. Eigen

`Eigen` is a C++ template (header-only) library for linear algebra. No installation is needed. Just download and decompress the source code:

```bash
$ wget -c https://gitlab.com/libeigen/eigen/-/archive/master/eigen-master.tar.gz -O - | tar -xz
``` 

This will create `eigen-master` under the current directory. It can be placed anywhere. To use this `Eigen`, 

<strong>NOTE</strong>: AxiSEM3D requires Eigen 3.3.9 or above, but the current stable version is 3.3.7 (up to July 22, 2020). Therefore, the above step is essential even you have installed `Eigen` using a package manager such as `conda` or `pip`. 

### 2. Boost
`Boost` provides free peer-reviewed portable C++ source libraries. `AxiSEM3D` uses 




## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTUxMTM0OTA2NCw0NjY4NzA2ODIsLTE2ND
cwNzg5MDksLTEzODM3NzAyMDYsLTE3NDkwNTg1MDUsMTM3MTg4
ODU4LC0zMzI3OTQ4NjcsLTE3Mzc1ODUxOTUsLTUyODkzNTk2MS
wxMTA3MDY4NjYwLC0yMTAwNDcxNjQ3LC0yMTYzMjEyMzgsMjIz
MDAyNzg1XX0=
-->