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
Eigen | linear algebra | 3.3.9 | The current stable release 3.3.7 (up to July 2020) is insufficient.
Boost | C++ helpers | 1.7.1 | `AxiSEM3D` only uses some of its header-only modules.
FFTW | fast Fourier transform | 3.3.8 | Both single- and double-precision builds are required.
Metis | mesh partitioning | 5.1.0 | Both 32- and 64-bit builds are accepted.
NetCDF | intensive structured IO | 4.7.1 | Parallel NetCDF is supported but not mandatory.





### 1. Installing dependencies locally
#### 1.1. Dependency

#### 1.1. Eigen

`Eigen` is header-only. Just download and extract the source code:

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
eyJoaXN0b3J5IjpbLTYzMzc3Njk2NCwtMTI3OTM1NDkxNCwxMj
E2MTk3MTQ1LC0xMzI3MDI2MjUwLC0xMzgxOTc0MzY4LDQ2Njg3
MDY4MiwtMTY0NzA3ODkwOSwtMTM4Mzc3MDIwNiwtMTc0OTA1OD
UwNSwxMzcxODg4NTgsLTMzMjc5NDg2NywtMTczNzU4NTE5NSwt
NTI4OTM1OTYxLDExMDcwNjg2NjAsLTIxMDA0NzE2NDcsLTIxNj
MyMTIzOCwyMjMwMDI3ODVdfQ==
-->