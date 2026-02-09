
# Element Output 

Rather than obtaining seismograms for specific stations, you may want to save seismic wavefields on more general portions of the mesh. AxiSEM3D allows you to do this with the element output, through which you can save wavefields on given elements through the mesh, specified by ranges. Similarly to station groups seen earlier, you can pass a list of element groups, e.g. defining azimuthal slices, or cuts at a given depth. 

<span style="color: red;"><b>*** TO DO: Fix link</b> </span>

Let us look at a real input file that you can find in [this example](https://github.com/kuangdai/AxiSEM-3D/blob/master/examples/03_Cartesian_SEG_EAGE_salt_5Hz/input1D/inparam.output.yaml) of a 5s run in a salt model on a cartesian mesh, and go through each option. 

You can see that in this example, `list_of_elements_group` contains two entries indexed by dashes: 

`- orthogonal_azimuthal_slices` and `- Fourier_coefficients_ocean_floor`. Each of those entries have the same parameters that need setting, and we will now go through them in detail. 

---

`elements`

You need to specify a horizontal and a vertical range which contains the elements. Note that an element will be considered within the range if its centre is.

-`horizontal_range`: The horizontal range, indexed by theta for a spherical mesh, and by s for a cartesian mesh.
-`vertical_range`: The vertical range, indexed by r for a spherical mesh, and by z for a cartesian mesh.

---

`inplane`

These are options relating to inplane properties of the elements, e.g. edges and GLL points.

-`edge_dimension`: Choose from "HORIZONTAL", "VERTICAL", or "BOTH". This specifies whether to record only on one outer edge of the element, or on both (i.e. the whole element).

-`edge_position`: ??

-`GLL_points_one_edge`: Allows you to only select a subset of GLL points to dump, and thus achieve inplane downsampling. For instance if you chose a simulation with a polynomial order npol = 4 for the spectral elements, you can pass [2] to only dump the central GLL point, [0,4] to only dump the vertices, [0,2,4] for the vertices, the edge centres, and the element centre; or pass "FULL" to dump everything.

---

`azimuthal`

So far we have only been specifying inplane properties, i.e. ranges and downsampling. Once this is done, you have the choice to dump all the Fourier coefficients for those elements, which allows you to create any azimuth you want; or you can only dump given azimuthal slices that you're interested in. This is accomplished by passing a list of azimuths, or a list of latitudes and longitudes. Note that if both lists are empty, you are de facto dumping all the Fourier coefficients.

-`phi_list`: A list of azimuths (in radians) for which the slices will be recorded.

-`lat_lon_list`: ???.  A list of latitudes and longitudes for the recorded slices. Azimuths computed from this list will be appended to `phi_list`.

---

`wavefields`

These parameters are the same as for the "Stations output" section.

---

`temporal`

These parameters are the same as for the "Stations output" section.

---

`file_options`

These parameters are the same as for the "Stations output" section, with the exception that the only available file format is NetCDF.

 

