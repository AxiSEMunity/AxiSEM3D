# inparam.model.yaml

## 1D model 
Parameters for 1D model (the mesh)

---

**Section:** model1D

**parameter name**:  exodus_mesh  
**what**: Exodus mesh file created by salvus mesher span  
**type**: filename.  
**default:**   
    
---

**Section:** geodesy  
Parameters for geodesy.  
**parameter name**:   lat_lon_north_pole_mesh. 
**what**:  geographic location of the north pole in the mesh.  
 **type**:  array of double / SOURCE.  
 **default:**   
Note:
1. This reference location enables the usage of geographic
coordinates for locating sources, receivers and 3D models,
compatible with Cartesian meshes 
2. array of double: [latitude, longitude]. 
3. SOURCE: determined this location by the FIRST source presented in list_of_sources in `inparam.source.yaml`; always use SOURCE for a single-source simulation.    


**parameter name**:   flattening_on_surface: WGS84. 
**what**:  flattening on the surface.  
**type**:  string / double.  
**only**: SPHERE / WGS84 / GRS80 / SPECFEM3D_GLOBE / value.  
**default:**     
Note: 
1) ellipticity is ignored for a Cartesian mesh
2) 0 for a perfect sphere; ~0.0033 for the Earth
3) ellipticity will be used in the transformation between
the geographic and the geocentric co-latitudes;
see eq. (14.32) in Dahlen & Tromp, 1998
4) to actually deform the entire mesh, add 3D model
"Ellipticity" to list_of_3D_models
    
---

**Section:** absorbing boundary  
Parameters for absorbing boundary condition

**parameter name**:  boundaries  
**what**:  model boundaries regarded as absorbing boundaries  
**type**:  array of string   
**only**: a subset of [RIGHT, BOTTOM, TOP].  
**default:**   
Note: 
1) an AxiSEM3D mesh may contain four outer boundaries:
left (axial), right, bottom and top (surface); the right,
bottom and top ones can be absorbing boundaries (the left
or axial one is non-physical)
2) use [] to disable absorbing boundary condition
  stress-free
3) the most common case in seismology is [RIGHT, BOTTOM]
  
**parameter name**:  enable_Clayton_Enquist
**what**:  enable the Clayton-Enquist approach.  
**type**: bool  
**default:**   
Note: the simplest linear approach by Clayton & Engquist (1977)


**Subsection:** Kosloff_Kosloff  
 The sponge approach by Kosloff & Kosloff (1986)
<div style="margin-left:48px;">  
<b>parameter name</b>: enable <br>
<b>what</b>: enable the Kosloff-Kosloff approach.  <br>
<b>type</b>: bool  <br>
<b>default:</b>  <br> 
Note: Clayton-Enquist and Kosloff-Kosloff can be used together, but one of them has to be enabled at least. 
<br><br> 
 
<b>parameter name</b>: relative_spans<br>
<b>what</b>: relative spans of the sponge layers<br>
<b>type</b>: array of double<br>
<b>default:</b>   <br> 
Note: 
</div>
<div style="margin-left:60px;">
1) must be presented in the same order as absorbing_boundaries<br>
2) to use Kosloff-Kosloff, the mesh should be a little largerthan the required computational domain; for example, givena required domain spans from 0 to 100 km in depth, one cangenerate a mesh from 0 to 110 km and set the relative span to 0.05, so the thickness of the sponge layer at the mesh bottom will be determined as 110 * 0.05 = 5.5 km, leaving an unaffected depth range from 0 to 104.5 km for normal wave propagation and analysis<br>
3) allowed range: .01 ~ 0.25<br>
</div>
<br>
<div style="margin-left:48px;">      
<b>parameter name</b>: gamma_expr_solid<br>
<b>what</b>: expression of γ-factor in solid domain<br>
<b>type</b>: math expression<br>
<b>default:</b>  <br> 
Note: </span><br>
</div>
<div style="margin-left:60px;">
1) γ-factor represents the absorbing strength at a point<br>
2) allowed arguments include (case sensitive):<br>
</div>
<div style="margin-left:70px;">
- VP, VS: P- and S- wave velocities at the point<br>
- RHO   : density at the po<br>
</div>
<div style="margin-left:60px;">
* VP, VS and RHO are the 1D values in the Exodus mesh<br>
3) this expression will be further multiplied by a pattern
function that equals to 1 on the outermost edge of the
sponge layer (i.e., on the mesh boundary) and gradually
decreases to 0 on the the innermost edge; such a decreasing
pattern is automatically handled by the solvers<br>
4) the default is an empirical expression from
Haindl et al., 2020<br>
</div>
<div style="margin-left:48px;">
<br>       
<b>parameter name</b>:   gamma_expr_fluid<<br>
<b>what</b>:  expression of γ-factor in fluid domain<br>
<b>type</b>:math expression<br>
<b>default:</b>  <br> 
Note: same as gamma_expr_solid but without VS dependency
</div>

---

**Section:** attenuation 

**parameter name**:  attenuation  
**what**:  attenuation mode  
**type**:  string. 
**only**:  NONE / FULL / CG4  
**default:**   
Note: 
1) NONE: turn off attenuation
2) FULL: compute attenuation on all GLL points
3) CG4:  compute attenuation on 4 GLL points per element;
CG4 is mostly as accurate as FULL but more efficient
than FULL, see van Driel & Nissen-​Meyer, 2014;
CG4 requires set(NPOL 4) in CMakeLists.txt;



---

**Section:** 3D models  

**parameter name**:  list_of_3D_models  
**what**:  list of 3D models.  
**type**:  array of objects.  
**default:**   
Note: 
1) the order in this list can affect the final 3D model
2) use [] if no 3D model presents

<span style="margin-left:48px;">**key**: arbitrary names</span><br>

<div style="margin-left:60px;">
<b>parameter name</b>:  activated<br>
<b>what</b>: activate this model<br>
<b>type</b>:  bool<br>
<b>default:</b>  <br> 


<b>parameter name</b>:  class_name<br>
<b>what</b>: class name<br>
<b>type</b>: string<br>
<b>default:</b>  <br> 
Note: current built-in classes include:
</div>  
<div style="margin-left:70px;">
 - StructuredGridV3D: volumetric 3D model on a structured grid<br>
 - StructuredGridG3D: geometric 3D model on a structured grid<br>
 - StructuredGridO3D: ocean-load 3D model on a structured grid<br>
 - Ellipticity: deform the mesh with global ellipticity<br>
</div>  
<br>

<div style="margin-left:60px;">
<b>parameter name</b>:  nc_data_file<br>
<b>what</b>: NetCDF data file for SStructuredGridV3<br>
<b>type</b>: filename<br>
<b>default:</b>  <br> 
<br>

<b>Subsection:</b> coordinates:<br>
Parameters for grid coordinates<br>
</div> 
<div style="margin-left:70px;">
<b>parameter name</b>: horizontal<br>
<b>what</b>:type of horizontal coordinates<br>
<b>type</b>: string<br>
<b>only</b>: DISTANCE_AZIMUTH / XY_CARTESIAN / LATITUDE_LONGITUDE<br>
<b>default:</b>  <br> 
<br>
<b>parameter name</b>: vertical  <br>
<b>what</b>: type of vertical coordinate<br>
<b>type</b>: string<br>
<b>only</b>: RADIUS / DEPTH<br>
<b>default:</b>  <br> 
<br>
<b>parameter name</b>: ellipticity: <br>
<b>what</b>: correct for ellipticity when locating the model<br>
<b>type</b>:  bool<br>
<b>default:</b>  <br> 
Note: used only when horizontal = LATITUDE_LONGITUDE<br>
<br>
<b>parameter name</b>:   depth_below_solid_surface<br>
<b>what</b>:  use solid surface as depth origin<br>
<b>type</b>: bool<br>
<b>default:</b>  <br> 
Note: used only when vertical = DEPTH<br>
<br>
<b>parameter name</b>: nc_variables<br>
<b>what</b>:  NetCDF variables for the coordinates<br>
<b>type</b>: array of string<br>
<b>default:</b>  <br> 
<br>
<b>parameter name</b>: data_rank<br>
<b>what</b>:   rank of the coordinates in data<br>
<b>type</b>: array of int<br>
<b>default:</b>  <br> 
<br>
<b>parameter name</b>: length_unit<br>
<b>what</b>:  length unit of the coordinates<br>
<b>type</b>: string / value<br>
<b>only</b>: km / m / number<br>
<b>default:</b>  <br> 
<br> 
<b>parameter name</b>:  angle_unit   <br>
<b>what</b>:  angle unit of the coordinates<br>
<b>type</b>:  string<br>
<b>only</b>: degree / radian<br>
<b>default:</b>  <br> 
<br>        
<b>parameter name</b>: undulated_geometry<br>
<b>what</b>: use undulated (otherwise reference) geometry to determine the vertical location<br>
<b>type</b>: bool<br>
<b>default:</b>  <br> 
Note: compatible only with vertical = RADIUS<br>
<br>
<b>parameter name</b>: whole_element_inplane<br>
<b>what</b>: check inplane model range for the whole element<br>
<b>type</b>: bool<br>
<b>default:</b>  <br> 
Note: <br>
</div>
<div style="margin-left:80px;">
1) if this parameter is set to true, the element center
will be used to determine whether an element is located
within the "inplane" model range<br>
2) if its center is in range, all its GLL points must be
in range, or an exception will occur; users can extend
the model range slightly to allow for numerical errors<br>
3) this parameter safely realizes inplane discontinuities<br>
<br>
</div>
<div style="margin-left:60px;">  
<b>Subsection:</b>   properties    <br>    
Parameters for properties<br>
</div>
<div style="margin-left:70px;">  
- VP:
</div>
<div style="margin-left:80px;">  
<b>parameter name</b>: nc_var <br> 
<b>what</b>: NetCDF variable <br> 
<b>default:</b>  <br> 
<br>
<b>parameter name</b>: factor        <br>         
<b>what</b>: factor or unit <br> 
<b>default:</b>  <br> 
<br> 
<b>parameter name</b>:                
<b>what</b>: reference kind<br>  
<b>only</b>: ABS / REF1D / REF3D / REF_PERTURB<br>  
<b>default:</b>  <br> 
Note: For any property X:<br>  
</div>
<div style="margin-left:90px;"> 
1) ABS: absolute value
 X_3D = value_in_file<br>  
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D<br>  
3) REF3D: perturbation w.r.t. the current 3D model
 X_3D = (1 + value_in_file) * X_3D<br>  
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D
reference_kind: REF1D<br>  
</div>
<br>  
<div style="margin-left:70px;">  
- VS:
</div>
<div style="margin-left:80px;"> 
<b>parameter name</b>:  nc_var<br> 
<b>what</b>: NetCDF variable<br> 
<b>default:</b>  <br> 
<br>   
<b>parameter name</b>:   factor<br> 
<b>what</b>:  factor or unit <br> 
<b>default:</b>  <br> 
<br>    
<b>parameter name</b>:   reference_kind   <br> 
<b>what</b>: reference kind <br> 
<b>only</b>: ABS / REF1D / REF3D / REF_PERTURB <br> 
<b>default:</b>  <br> 
Note: For any property X: <br> 
</div>
<div style="margin-left:90px;">  
1) ABS: absolute value X_3D = value_in_file <br> 
2) REF1D: perturbation w.r.t. the 1D reference modelvX_3D = (1 + value_in_file) * X_1D <br> 
3) REF3D: perturbation w.r.t. the current 3D model X_3D = (1 + value_in_file) * X_3D <br> 
4) REF_PERTURB => perturbation w.r.t. the current
p X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D  <br> 
</div>      
<br> 
<div style="margin-left:60px;">   
<b>parameter name</b>:  store_grid_only_on_leaders<br> 
<b>what</b>: store grid data only on the leader processors<br> 
<b>type</b>: bool<br> 
<b>default:</b>  <br> 
Note: turn this on if the model is large; set mpi:nproc_per_group
in inparam.advanced.yaml to the number of processors per
node to minimize memory usage
</div>