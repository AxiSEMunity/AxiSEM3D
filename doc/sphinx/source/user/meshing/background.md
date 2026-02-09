# Background Files (.bm)

As mentioned earlier, creating a mesh requires a model file. These are saved in a plain text format which we denote .bm (for background model). At the simplest level, the .bm file just needs to contain the variation of seismic parameters (density, v_p and v_s) with depth, or equivalently, radius. You can edit .bm files using any old text editor, but you might find Visual Studio's functionality useful for this, if you already have it installed. 

The first five lines of the .bm file contain header information (ignore bullet points): 

* NAME - states what the model is called. You can set this however you like 
* ANELASTIC - states whether the model contains attenuation in the .bm file [T/F or True/False]
* ANISOTROPIC - states whether your model contains seismic anisotropy in the .bm file [T/F or True/False]
* UNITS - states the units of parameters. [Metres\/kilometres]
* COLUMNS - a list of column headers. [At the very least, should include radius or depth, vp, vs, and rho (density)]

Note that your file must be consistent: if you create or use a .bm file with `ANISOTROPIC` and `ANELASTIC` set to `T` (i.e. your model is both attenuating and non-isotropic), you need to include the correct column headers and correct corresponding values. This means that an model with `ANISOTROPIC = T` must have values for vpv and vsv as well as vph and vsh. Similarly, if you set `ANELASTIC' to 'T', make sure to include 'QKa' (Q_Kappa), 'QMu', and 'eta' in your column headers and seismic profiles. If you do not do this, and try and generate a mesh for which the headers and columns are inconsistent, the mesher will throw an error. 

Here's an example of a martian .bm file which is both attenuating and anisotropic: 

* NAME SEIS_shallow_layer
* ANELASTIC T
* ANISOTROPIC T
* UNITS m

|COLUMNS |radius |rho| vpv |vsv |qka| qmu |vph |vsh| eta|
|--- |--- |---| --- |--- |---| --- |--- |---| ---|
|| 0. | 6776.00 |5635.00 |0.00 |57822.0 | 143.0 |5935.00  | 0.00 |1.00000|
|| 19303.| 6776.00| 5635.00 |0.00| 57822.0|  143.0 |5935.00  | 0.00| 1.00000|
||38607.| 6776.00 |5635.00|0.00 |57822.0 | 143.0 |5935.00  | 0.00 |1.00000|

If you want to include fluid regions, you just do this by setting the v_s values to be zero. You can have fluid layers at the surface if you are running simulations with, for example, an ocean. This example has only three rows (from 0 to 38km, i.e. all within the martian core), but you can have many thousands of rows, if you so wish. The mesher will interpolate seismic parameters between your given values, and design an appropriate mesh. 

Note that in a .bm file, discontinuities are designated by a two consecutive lines with the same radius/depth value, but different entries in some or all of the other columns. This will place an element boundary at exactly the discontinuity location. You may, for certain purposes, wish to force an element boundary into a certain radius in the mesh. To do this, just repeat a line with all entries the same. 

There are other formats that you can save your .bm file in, but we'll stick to this one here because it's the easiest to understand. Note, on some systems, you may get mesher errors if you do not indent all the rows containing numerical values by two leading tab spaces.  
