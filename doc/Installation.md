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

The solver is developed on top of several modern numerical packages including `fftw`


## Tools for pre- and post-processing




[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTY4Mjk1ODQ3NSwtNTI4OTM1OTYxLDExMD
cwNjg2NjAsLTIxMDA0NzE2NDcsLTIxNjMyMTIzOCwyMjMwMDI3
ODVdfQ==
-->