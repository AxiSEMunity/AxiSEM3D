
# Using Volumetric Models
<span style="color: red;"><b>*** TO DO: This is a mashup of 2 different sources so it needs editing. It is also very long and probably should be sectioned into different files.</b> </span>

The 3D medium simulated by AxiSEM3D has the shape of a body of revolution with the cross section prescribed by the 2D mesh. Using a local-scale simulation as an example: if the 3D input model is given in Cartesian coordinates, i.e. by nature in shape of a cuboid, AxiSEM3D only reads the subset of this model which fits into the cylinder which is being simulated. 3D models which are smaller than the simulated volume may also be used. In this case regions not covered by the 3D model are filled in with the background model given in the 2D mesh.

The medium properties which can be set using volumetric models are:
* homogeneous: "VS", "VP", "RHO"
* transversely isotropic: "VSV", "VSH", "VPV", "VPH", "ETA"
* fully anisotropic: (this is the stiffness tensor in voigt notation; all components must be provided)

    "C11", "C12", "C13", "C14", "C15", "C16", 
    
    "C22", "C23", "C24", "C25",  "C26", 
    
    "C33", "C34", "C35",  "C36", 
    
    "C44", "C45",  "C46", 
    
    "C55",  "C56", 
    
    "C66" 

Once you have created a new 3D model and saved it as a netcdf file, you
need to upload it to the input folder. Then, there are two input files
that you need to alter in order to get the model to be read in and
represented properly: `nparam.model.yaml` and `inparam.nr.yaml`.

## inparam.model.yaml

The portion of the `inparam.model.yaml` file that needs to change from
the 1D case is the `3D Models` section. Under the header `list_of_3D models`, you can add as many as you want. Let us start with `etopo1`,
which you can put as the `key` in the first indented line. If you are
doing checks and want to run the code without the model being read in,
but don’t want to rewrite the whole input file, just toggle `activated`
to `false`.

`class_name` for all 3D geometric models will be `StructuredGrid3D`.

`class_name` for a 3D volumetric model on a structured grid will be
`StructuredGridV3D`.

`nc_data_file` is also just the name of the 3D model file that you have
created and saved in the /input folder.

For coordinates, this project will probably use `LATITUDE_LONGITUDE` and
either `DEPTH` or `RADIUS`. Just make sure you get the latter correct:
`etopo1` is defined as an elevation with respect to some mean sea level
whereas the Moho file we are using is defined as a depth with respect to
sea level.

Again, you can use the `ellipticity` correction if you want to but it is
unlikely to make an enormous difference.

If you’ve used the provided python scripts, you should not need to
change the `data_rank`, `length_unit`, `angle_unit`, or `factor` flags
from those we have given (*\[0,1\], m, degree* and *1.0* respectively).

The things that you will need to change are as follows. Close and
careful attention is needed when changing these, as mistakes can be very
challenging to debug. The order that we go through these is not the
order that they appear in the input file, but we think it is more
logical.

- `nc_variables` and `nc_var`: these are the variable names used in the
  netcdf file. If you have kept the defaults, the first one should be
  *\[latitude, longitude\]* and the second one should be *depth*, moho
  and elevation for etop01 and *C*<sub>11</sub>, *C*<sub>12</sub>... for
  seismic anisotropy. Note that these are just the labels that we used
  in our netcdf file, they can be called anything and just because you
  set ‘depth’ or ‘radius’ as the name does not mean that the code will
  treat it as such unless you also set the flags below.

- `vertical` should be `DEPTH` for the Moho and seismic anisotropy, and
  RADIUS for topography. Make doubly sure this is consistent, if you
  have files which you think are the wrong way round and do not want to
  remake them, you can always swap `factor` from *1.0* to *-1.0*.

- `undulation_range: interface` should be the depth or radius of the
  interface you are undulating, as defined in the background 1D model.
  If you are using depth as your variable (as we used for the Moho),
  then this should be, for example, -*24400.0*, as the Moho is 24.4 km
  below the surface in PREM. If you are using radius as your variable,
  as we do for the surface topography, this should be *6371000.0*.
  Remember that you will define the undulations as elevations above the
  boundary in consideration. So, at a subduction trench where the Moho
  is at 7 km below the surface, you would define the gridpoint for the
  Moho’s depth as +17400 (-7000 + 24400 = 17400 m). Similarly, in Tibet
  where the Moho is at 80 km depth, the gridpoint would be -55600
  (-80000 + 24400 = -55600 m).

* `depth_below_solid_surface`: This option is only relevant if there is an ocean or other fluid medium at the free surface and a solid medium below, and if depth is used as a vertical coordinate. In this case, if set to true, depth will be measured as depth below the ocean floor rather than as depth below the ocean surface.

    Lots of models in
  seismology (e.g., PREM) do not include the ocean, so if you add one on
  top, you need to either add one right on top (making the total radius
  6371 km), or remove 3km of rock first (to keep the radius as 6368 km).
  This flag allows you to define whether your interface coordinates are
  with respect to the seafloor (`true`) or the ocean surface (`false`).
  You might also want to think about whether, when you read in a 3D
  model for subsurface features, the depths should be defined with
  respect to the solid surface in the 1D case, or the 3D case – as the
  position of the surface may have moved considerably in between these.

* `data_rank`: This option allows for flexible handling of row-major and column-major input data. AxiSEM3D is written in C++ and hence uses column-major data storage by default. If the input data is also in column-major format (e.g. default python output) the data rank option should be set to *[0,1,2]*. If the data is in row-major format (e.g., default MATLAB and Fortran output) the order needs to be reversed, to *[2,1,0]*. This option also allows the user to easily swap axes if they need to debug their model setup.

* `undulated_geometry`: This option is only relevant if geometric models (including ellipticity) and volumetric models are used at the same time, and it determines the order in which the models are read. As explained in the next section, geometric models effectively deform the mesh. If this option is set to `false`, the volumetric model is read in first and then deformed together with the mesh. If the option if set to `true`, the mesh is deformed first and the volumetric model is read afterwards, such that it is unaffected by the geometric model.

- `undulation_range: min_max`
   tells the code in which
  range it is permissible to stretch or shrink elements to accommodate
  your deformed boundary. For example, if you set the Moho to be at
  24 km and then undulate it, elements on either side need to deform to
  accommodate the undulated boundary. `min` and `max` define the lower
  and upper edges of the deformation zone, which should be at least
  twice as large as the undulation range. For example, in Crust 1.0 the
  Moho depth varies between 75 km and 7 km, so doubling the deepest
  depth and halving the shallowest one in this range would give you
  *min_max = \[3500.0, 150000.0\]*. Again, just make sure that the signs
  are consistent and correct: 3500 m is actually a shallower depth
  (i.e., the minimum deformed depth), rather than the maximum! Either
  way the interface needs to be contained within the bounds you have
  set, else the undulation will not work!

* `whole_element_inplane`: This option is used to make 3D models adhere to the mesh geometry. If it is set to `false`, each point inside an element reads the volumetric model individually. If set to `true`, it is ensured that the model boundaries are placed at the edge of elements. Use this option if a 3D model borders on a solid-fluid boundary or a sharp velocity contrast.

For volumetric models, `nc_variables*`should be set in the order *\[lat,
lon, depth\]* or *\[lat, lon, radius\]*. The names should match the
names in the netcdf file.

Data rank should map to the order of the variables in the netcdf file.
For example, if the netcdf file is in the order radius, longitude,
latitude, then the data rank would be set as follows: *data_rank: \[2,
1, 0\]*.

> **Note on Efficiency.**  
>Before starting the time-loop, AxiSEM3D attempts to locally reduce 3D input models to 1D as well as to decrease the level of anisotropy. Hence, the performance of the simulation is automatically optimised and does not depend on the structure of the input data.
>
>However, large volumetric models are a bottleneck in terms of working memory, since the entire model has to be read once before being distributed to the parallel nodes. If the model can be divided up into separate sections with the meshed background model in between this is highly recommended.

## inparam.nr.yaml

### Fourier series introduction

This input file, which we did not touch upon earlier, is used to set the
azimuthal (i.e. variation in longitude) parameters of the simulation.
You will recall that AxiSEM3D works by taking a 2D mesh and representing
the third dimension through a pseudospectral (Fourier) expansion. The
order of this expansion determines the complexity of the entire
simulation in the azimuthal direction. Note that this does not just
affect the degree to which 3D models are accurately represented, lower
Fourier order expansions also reduce the azimuthal complexity of the
seismic wavefield that you are going to record.

To understand this another way, the Fourier expansion in AxiSEM3D is
akin to the representation of any normal mathematical function in terms
of sines and cosines. The more terms you include, the more accurate the
function’s approximation is. For an arbitrary function, an infinite
number of terms is needed for an exact solution - but we can make do
with a smaller (finite) number of terms *Nu,* *Nr* if we are happy for
the solution to only be accurate to within a given degree. *what is the
conversion between these two? Not sure it is consistently defined...*

How much smaller than infinity the number of terms/Fourier order *Nr*
can be depends on a number of different factors:

1.  The complexity of the 3D model that you add to AxiSEM3D. If you have
    a very long-wavelength model with relatively smooth variations, like
    SEMUCB, you need far fewer terms than if you have something that is
    strongly varying and densely sampled, like ETOP01. For a 1D model,
    you never need more than *Nr = 5*.

2.  The complexity of the wavefield that you want to represent: in
    general, you will only add in complex 3D models if you want to see
    what impact they have on synthetic seismograms. However, if you
    wanted for some reason to use a complex 3D model but were not
    interested in the fine detail of the wavefield, you could get away
    with a lower Fourier order.

3.  The complexity of the source: as noted before, if you are only using
    a 1D model you can save a small amount of time by reducing the
    Fourier order needed to represent the source term. A second rank
    **moment tensor** (i.e. something from CMT) requires `Nr = 5`, a **point
    force** requires `Nr = 3`, an **isotropic/monopolar** pressure term
    requires `Nr = 1`. For most 3D models the value of *Nr* needed to
    capture their variations will be far larger than 5, so there is no
    saving to be made even if you only use a point force or pressure
    term.

Remember that the overall resolution of the simulation (i.e., the
smallest scale that you can resolve meaningfully) is still limited by
the mesh resolution, so increasing the Fourier order too much becomes
pointless as you are not sensitive to these small spatial scales anyway.
For this reason you should leave `bound_Nr_by_inplane = true`.

### Choosing input parameters

You will probably be wondering by now, ‘how small a Fourier order can I
get away with?’ - the answer is of course that this depends on how
accurate a simulation you want, and what much computational resources
you have available - a higher Fourier order increases both runtime and
memory demands. There are four options that you can choose in AxiSEM3D
for how the Fourier order is calculated, which may help with this
question. These are, `type_Nr =` :

- `CONSTANT`: here, the same *Nr* is used throughout the mesh. As
  mentioned previously, for a 1D model with a moment tensor 5 is both
  necessary and sufficient.

- `ANALYTICAL`: this is like constant, but with slightly finer
  user-control, in that you can choose how Nr varies through the
  interior of the planet. As an example, if you wish to investigate the
  effects of crustal structure on surface waves, you might want a high
  expansion in the mantle, but not in the core. This is thus a finer
  level of control than in *CONSTANT*, which can save you both time and
  memory.

- `STRUCTURED` allows you to specify even finer degrees of control, by
  setting the fourier order on individual elements. You can thus save
  even more time and resource with respect to *ANALYTICAL*, but a more
  detailed knowledge of how to create the appropriate grid of points is
  needed.

- `POINTWISE` allows you to go even further by having the code try and
  optimise the variation of *Nr* through the mesh. See below for more
  details.

Within `ANALYTICAL`, you can specify the variation of *Nr* with depth
that you wish to use within the `depth_dependent_AxiSEM3D_default`
option. Using the options `control_depths` and `Nr_at_control_depths`
you can make a pointwise specification of how Nr varies, remembering to
give your depths in metres.

For example, `control_depths: \[0., 6371.\]` and 
`Nr_at_control_depths:\[100, 100\]` would set `Nr = 100` across the whole of the Earth’s
radius. Changing `Nr_at_control_depths` to *\[100, 1000\]* would set
*Nr* to *100* at the surface, and *1000* at the planet’s centre, with
linear interpolation in between. You can add as many control points as
you want into these arrays. For example, if you wanted to look at
high-resolution crustal structure, you might do something like:

`control_depths: \[0, 50e3, 200e3, 6371e3\]` and 
`Nr_at_control_depths:\[1000,1000, 100, 10\]` which would give you 
`Nr = 1000` in the top
*50* km, decreasing to *100* in a linear fashion by *200* km depth, and
then linearly again from 100 at *200* km to 10 the centre.

### Wavefield scanning

Wavefield scanning is an advanced feature of AxiSEM3D. It involves the
code trying to ‘learn’ the most optimal *Nr* profile it can, using past
runs as a guide. In a way, this lets the code do all the optimisation
described in the steps above itself, such that the runtime and memory
requirements are minimised.

This method works by analysing the frequency spectrum at all the
interpolation points in the mesh. For a given point, if the energy above
a given frequency is below some threshold value (which you can set), the
code determines that frequencies above this point are not contributing
to the overall solution, and the Fourier series can thus be truncated.
In subsequent runs, the value of *Nr* at these points can therefore be
reduced as they are determined to have a reduced azimuthal wavefield
complexity.

If you want to try this, you can set `enable_scanning` to `true` and
`type_Nr = POINTWISE`. Follow the instructions in the input file to make
sure that you match the filenames and copy the relevant files across.