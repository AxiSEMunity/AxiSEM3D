# Time Axis

This sets the length of your simulation, the timestep (if you want to set it manually), the Courant number, and the integrator type. The descriptions in the file are pretty comprehensive, a casual user will only need to change the record length. Note that the simulation will be slightly longer than record_length, due to the 'ramp-up' of the source time function at the start. You can shift or crop the end result if this is a problem, or use a delta function source (described below) which will minimise the amount of extra added time. 

# Sources 

Here, you can set the number of sources you want to use. Using no source is only recommended for efficiency/scaling tests, as without a source, nothing happens, and your wiggly lines become flat. The input files are relatively self-explanatory, and most users will be concerned only with the choice of source location and source-time function shape/duration. The input parameters are divided into three categories: 

* Location
* Mechanism 
* Source-time Function

The main difference from the previous releases of AxiSEM3D and AxiSEM is that here, we do not need to perform a coordinate rotation such that the source lies on the axis. This means that now you may have more than one source, and these sources do not need to lie on a straight line through the centre of the planet. This means that `Location` is the category with the biggest differences to previous releases, especially if you have multiple sources. 

#  Mechanism and Source-time function

Your source can be either a zeroth, first, or second rank tensor: corresponding to an isotropic pressure source, a point force vector, or a full moment tensor. Bear in mind that the moment tensor follows the `CMTSOLUTION` convention, which uses the non-SI unit of dyne-centimetres. 

If you are used to using other codes, just make sure you check the order of the components in the moment tensor. These, and the coordinate system, are clearly identified in the input file documentation. 

You can also choose whether to use a Gaussian source time function, or a user-defined one that you want to read in using `StreamSTF` or `netCDF_STF`. 

If you want to use a delta function, set the half-duration to zero (the code will define the half duration as being 5 x dt, which is much less than the mesh period). For some purposes (e.g. wavefield visualisation), you need to set the source-time function parameters before you do the simulation. For 'ordinary' seismograms, it's fine to use a delta-function source and deal with the rest in post-processing through a convolution. This is probably more efficient for lots of purposes, as changing the source-time function will not mean re-doing the entire simulation. 

# Location and Multiple sources 


You have two options for how to choose your coordinate system. Note that these are set in the `inparam.model` file, not `inparam.source` - but you need to be consistent across the two files. The coordinate system in the model file is determined by your choice of axis. 

Your options are: 

```
lat_lon_north_pole_mesh = source
```
or 
```
lat_lon_north_pole_mesh = [latitude, longitude]
```

The first case means that the FIRST source position is given in latitude-longitude, and then this will become the axis used in simulations (again, for details of how this works see papers, but in short it sets the geographic location of the symmetric axis). If you're only using one source this is a good default and you should be careful about changing the coordinate system and unnecessary coordinate rotations are a recipe for making mistakes. 

The second case means that you give the geographic location of the symmetric axis. If you do this, you can give the source coordinates in either latitude-longitude or distance-azimuth (note that the depth parameter does not change either way). The sources are all then off-axis, unless you deliberately place them on the axis. 

The second case is advantageous if you have a situation where you don't have a series of sources all along a radial line. This might include distributed sources, such as finite fault approximations or adjoint simulations. 

