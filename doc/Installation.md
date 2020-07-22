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

This will create `eigen-master` under the current directory. Next, in `SOLVER/CMakeLists.txt` 



## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTk0ODQ5OTY4OCwtMzMyNzk0ODY3LC0xNz
M3NTg1MTk1LC01Mjg5MzU5NjEsMTEwNzA2ODY2MCwtMjEwMDQ3
MTY0NywtMjE2MzIxMjM4LDIyMzAwMjc4NV19
-->