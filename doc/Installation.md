[<< Back to repository](https://github.com/kuangdai/AxiSEM-3D)


# Installation

The installation of AxiSEM3D includes three parts: the mesher, the solver and several tools (mostly python libraries) for pre- and post-processing. 


## Mesher
[SalvusMeshLite](https://gitlab.com/Salvus/SalvusMeshLite) is the mesher for AxiSEM3D. Its installation is trivial with `pip`: 

```bash
$ pip install https://gitlab.com/Salvus/SalvusMeshLite/-/archive/master/SalvusMeshLite-master.zip
```

Verify the installation by

```bash
$ python -m salvus_mesh_lite.interface AxiSEM -h
```

This will display all the arguments you can pass to the mesher. 



## Solver

## Tools for pre- and post-processing
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTIxMDA0NzE2NDcsLTIxNjMyMTIzOCwyMj
MwMDI3ODVdfQ==
-->