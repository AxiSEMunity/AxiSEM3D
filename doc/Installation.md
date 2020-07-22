[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation

The installation of AxiSEM3D includes three parts: the mesher, the solver and several tools (mostly python libraries) for pre- and post-processing. 


## Mesher
[SalvusMeshLite](https://gitlab.com/Salvus/SalvusMeshLite) is the mesher for AxiSEM3D. Its installation is trivial with `pip`: 

```bash
$ pip install https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip
```

There is no need to install the mesher on an HPC cluster. The mesher runs on a single processor, but is efficient enough to generate large-scale meshes in 

Verify the installation by

```bash
$ python -m salvus_mesh_lite.interface AxiSEM -h
```

which will display all the arguments you can pass to the mesher. 


## Solver

## Tools for pre- and post-processing
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTg5MjQ1OTg2MiwtMjEwMDQ3MTY0NywtMj
E2MzIxMjM4LDIyMzAwMjc4NV19
-->