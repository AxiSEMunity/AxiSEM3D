[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation

The installation of AxiSEM3D includes three parts: the mesher, the solver and several tools (mostly python libraries) for pre- and post-processing. 


## Mesher

[SalvusMeshLite](https://gitlab.com/Salvus/SalvusMeshLite) is the mesher for AxiSEM3D. Its installation is trivial with `pip`: 

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

The AxiSEM3D solver is developed on top of several modern numerical packages including `Eigen`, `Boost`,  `FFTW`, `Metis` and `NetCDF`. 

### Eigen

Eigen is a C++ template (header-only) library for linear algebra. We only need to download and decompress the source code:

```bash
$ wget -c https://gitlab.com/libeigen/eigen/-/archive/master/eigen-master.tar.gz -O - | tar -xz
``` 

This will create `eigen-master` under the current directory. 

<strong>NOTE</strong>: AxiSEM3D requires Eigen 3.3.9 or above, but the current stable version is 3.3.7 (up to July 22, 2020). Therefore, the above step is essential even you have installed Eigen using package managers such anaconda 




## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTEwMzQ2MDMxNzEsMTM3MTg4ODU4LC0zMz
I3OTQ4NjcsLTE3Mzc1ODUxOTUsLTUyODkzNTk2MSwxMTA3MDY4
NjYwLC0yMTAwNDcxNjQ3LC0yMTYzMjEyMzgsMjIzMDAyNzg1XX
0=
-->