# 3D Models

AxiSEM3D differentiates between two types of 3D models: **volumetric** and **geometric**. The volumetric model is provided in the form of 3D grid data which contain medium properties of a part of or the entire modeled volume. The mesh geometry remains unchanged by this. Geometric models, on the other hand, can be used to impose undulations on the mesh. This facilitates the simulation of non-horizontal solid-fluid interfaces, and strong velocity contrasts can be modeled accurately.

Ellipticity is a hard coded special case of a geometric model. It is also possible to enable Ellipticity as part of the coordinate transform without deforming the mesh.

It is not possible to use a volumetric model to impose a solid medium onto a region which is fluid in the 2D mesh or vice versa. All fluid and/or solid regions currently need to be continuous. However, AxiSEM3D offers a 3D ocean load function which allows approximate modeling of discontinuous oceans and continents.

Sufficient azimuthal sampling is crucial for the proper modeling of 3D media. It is therefore highly recommended to read the Nr Field section before running a big simulation with 3D models.

---

```{toctree}
---
maxdepth: 2
---
geometric.md
volumetric.md
volumetric-usage.md
ocean.md
absorbing.md
debugging.md
```
