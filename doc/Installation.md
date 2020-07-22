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


### 1. Eigen

`Eigen` is header-only, so we just need to download the source code:

```bash
$ git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
```
This will create `eigen3_develop` under the current directory. To let AxiSEM3D find it, do

```bash
$ export EIGEN3_ROOT=$PWD/eigen3_develop
``` 

<strong>NOTE</strong>: `AxiSEM3D` requires `Eigen` 3.3.9 or above, but the current stable release is 3.3.7 (up to July 22, 2020). Therefore, the above steps are essential even you have installed `Eigen` before. 


### 2. Boost
`AxiSEM3D` only uses some of the header-only modules of `Boost`, so we just download the source code:

```bash
$ git clone https://gitlab.com/libeigen/eigen.git eigen3_develop
```




## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTM4NjA4Mjg1NywxMjE0MDIxMjIsLTE5Mz
I5MjQyNzYsLTYzMzc3Njk2NCwtMTI3OTM1NDkxNCwxMjE2MTk3
MTQ1LC0xMzI3MDI2MjUwLC0xMzgxOTc0MzY4LDQ2Njg3MDY4Mi
wtMTY0NzA3ODkwOSwtMTM4Mzc3MDIwNiwtMTc0OTA1ODUwNSwx
MzcxODg4NTgsLTMzMjc5NDg2NywtMTczNzU4NTE5NSwtNTI4OT
M1OTYxLDExMDcwNjg2NjAsLTIxMDA0NzE2NDcsLTIxNjMyMTIz
OCwyMjMwMDI3ODVdfQ==
-->