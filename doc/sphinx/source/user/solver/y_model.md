# inparam.model.yaml

This is probably the file that you will spend the most time editing.

In this file, you specify:

**1D model**

 The name of your exodus mesh file which should end in in “.e".

```
model1D:
    exodus_mesh: global_mesh__prem_ani__50s.e
```

---

**geodesy**

The location of the north pole in the mesh and flattening in the geodesy section. For this example, the pole is determined by the location of the FIRST source
presented in `list_of_sources` in `inparam.source.yaml` and uses the WGS84 geodetic model.

```
geodesy:
    lat_lon_north_pole_mesh: SOURCE
    flattening_on_surface: WGS84
```

---
**absorbing boundary**

Describes the absorbing boundary conditions. For this example, the absorbing boundaries are the right and bottom boundaries, the Clayton & Engquist (1977) and the Kosloff & Kosloff (1986) approach are enabled. The relative span of the sponge layer is specified along with expressions for the γ-factor in solid  and fluid domains

```
absorbing_boundary:
    boundaries: [RIGHT, BOTTOM]
    enable_Clayton_Enquist: true
    Kosloff_Kosloff:
        enable: true
        relative_spans: [.05, .05]
        gamma_expr_solid: 1.1 / T0 * (VS / VP)^2 * exp(-0.04 * SPAN / (VP * T0))
        gamma_expr_fluid: 0.88 / T0 * exp(-0.04 * SPAN / (VP * T0))
```
---

**attenuation**

Attenuation is computed on 4 GLL points per element.

```
attenuation: CG4
```
---


**3d models**

Since this is a 1D model, no 3D model is specified.
```
list_of_3D_models: []
```

For more details, see the section on 3D Models.