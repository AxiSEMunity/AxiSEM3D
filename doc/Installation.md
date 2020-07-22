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

The `AxiSEM3D` solver is developed on top of several modern numerical packages including `Eigen`, `Boost`,  `FFTW`, `Metis` and `NetCDF`:


Name|Usage|Version
--- | --- | ---
Eigen | linear algebra | **nicely**
1 | 2 | 3


### Installing dependencies locally



### 1. Eigen

`Eigen` is a C++ template (header-only) library for linear algebra. No installation is needed. Just download and decompress the source code:

```bash
$ wget -c https://gitlab.com/libeigen/eigen/-/archive/master/eigen-master.tar.gz -O - | tar -xz
``` 

This will create `eigen-master` under the current directory. It can be placed anywhere. To enable `AxiSEM3D` to find this `Eigen`, do

```bash
$ export EIGEN3_ROOT=$PWD/eigen-master
``` 

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the current stable version is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even you have installed `Eigen` using a package manager such as `conda` or `pip`. 


### 2. Boost
`Boost` provides free peer-reviewed portable C++ source libraries. `AxiSEM3D` uses header-only part of `Boost`, so we only need to download  




## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMzYzMDM3OTQxLC0xMzI3MDI2MjUwLC0xMz
gxOTc0MzY4LDQ2Njg3MDY4MiwtMTY0NzA3ODkwOSwtMTM4Mzc3
MDIwNiwtMTc0OTA1ODUwNSwxMzcxODg4NTgsLTMzMjc5NDg2Ny
wtMTczNzU4NTE5NSwtNTI4OTM1OTYxLDExMDcwNjg2NjAsLTIx
MDA0NzE2NDcsLTIxNjMyMTIzOCwyMjMwMDI3ODVdfQ==
-->