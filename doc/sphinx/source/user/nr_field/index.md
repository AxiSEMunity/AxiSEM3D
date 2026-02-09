# Nr Field
The Nr field is the heart piece of the AxiSEM3D method. It sets it apart from standard 3D methods, giving it versatility and the potential to save significant computational cost. Moreover, the complexity maps obtained from [wavefield learning](learning.md) can be used to analyse and describe the wavefield. The Nr field functionality may seem daunting at first, but it increases the control the user has over the cost-accuracy trade-off and it simplifies the progression from cheap test simulations to full simulations.

As explained in the meshing section, AxiSEM3D uses a 2D mesh despite being a fully 3D method. Sampling in the 3rd dimension is determined by the Nr field. Specifically, Nr is the number of evenly spaced sampling points in the azimuthal direction of the cylindrical coordinate system. The number can be different for each point on the 2D mesh (i.e. Nr(s,z)). To retain roughly constant sampling intervals, Nr should increase linearly with distance from the axis.

---

```{toctree}
---
maxdepth: 2
---
running.md
optimizing.md
learning.md
debugging.md
```
