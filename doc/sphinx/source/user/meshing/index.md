# Meshing

## Overview
AxiSEM3D makes use of the stand-alone mesher [SalvusMeshLite](https://gitlab.com/swp_ethz/public/SalvusMeshLite). 

See Installation instructions for [The Mesher](../../installation/mesher.md). 

## Usage

Please see the [SalvusMeshLite](https://gitlab.com/swp_ethz/public/SalvusMeshLite) repository for the latest instructions.

To create a self-explanatory `input.yaml` YAML file, run: 

```bash
$ python -m salvus_mesh_lite.interface AxiSEM --save_yaml=input.yaml
```

Edit it to your satisfaction and run:

```bash
$ python -m salvus_mesh_lite.interface --input_file=input.yaml
```

This creates and Exodus mesh file (*.e) required by the solver, AxiSEM3D.

Remember that all AxiSEM3D meshes are two-dimensional (i.e. flat) and use a one-dimensional seismic profile (i.e. there is no lateral variation in seismic parameters across the mesh) - 3D features are added in separately from the mesh.

You can visualise the meshes you create in python, but it is easier to use [Paraview](https://www.paraview.org) - which you can download for free. More on this later. 

The mesher is extremely efficient, and you should be able to produce fine meshes, even on a global scale, using a laptop or desktop computer. There is therefore no need to install the mesher on your HPC architecture; rather you can just copy across the finished mesh which will be saved as an `Exodus` (.e) file. 

## The two different types of AxiSEM3D meshes

Typing the line above should show you two available commands within the mesher which are relevant to AxiSEM3D: `AxiSEM` and `AxiSEMCartesian`. These are the two types of meshes that are compatible with AxiSEM3D: a global mesh (i.e. one that uses a curved geometry) and a cartesian one (which can have curvature on the surface, but is not designed for global-scale use). You can read the AxiSEM3D papers for more details on the differences between these meshes, but here we will focus on how to use them. 

---

```{toctree}
---
maxdepth: 1
---
background.md
gen_mesh.md
```
