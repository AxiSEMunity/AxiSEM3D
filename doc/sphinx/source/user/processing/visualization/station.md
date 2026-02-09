# Station Output

A standard reason to run a forward model is to obtain seismograms at given locations on the planet, e.g. corresponding to the locations of real seismic stations. In AxiSEM3D you can pass a list of so called station groups that will define the desired locations and outputs. A station group can for instance be the global seismic network containing all the corresponding stations, or the US array, or even a collection of stations bearing no connection to reality. 

<span style="color: red;"><b>*** TO DO: Fix link</b> </span>
Let us look at a real input file that you can find in [this example](https://github.com/kuangdai/AxiSEM-3D/tree/master/examples/01_spherical_Earth_PREM_50s/input/inparam.output.yaml) of a 50s run in PREM, and go through each option. 

You can see that in this example, `list_of_stations_group` contains two entries indexed by dashes: 

`- global_seismic_network_GSN` and `- USArray_transportable`. Each of those entries have the same parameters that need setting, and we will now go through them in detail. 

**`locations`**

- `station_file`

For this we need to specify a `station_file`, which is a text file where each row is a station, and the columns are "name", "network", "x1", "x2", "useless", "x3". "useless" is only present for reasons of compatibility with specfem, and the values in that column are discarded. "x1", "x2", and "x3" can take on different values, viz. "latitude", "distance", or "x" for x1, "longitude" or "azimuth" for x2, and "radius" or "depth" for x3. The meaning of each of these options is detailed in the next two points.
- `horizontal_x1_x2`:
Choose from "LATITUDE_LONGITUDE", "DISTANCE_AZIMUTH", or "XY_CARTESIAN". "LATITUDE_LONGITUDE" are the standard way of specifying stations in a geocentric system, with units of radians. "DISTANCE_AZIMUTH" is a way to specify stations in a source-centered coordinate system, where "distance" is in meters for cartesian meshes and radians for spherical meshes, and azimuth is in radians. "XY_CARTESIAN" is specifically for cartesian meshes, and is in meters.

- `vertical_x3`:
Choose from "RADIUS" or "DEPTH". "RADIUS" counts from the center of the Earth, and "Depth" from the surface. For both the unit is meters.

- `ellipticity`:
If `horizontal_x1_x2` is set to "LATITUDE_LONGITUDE", setting `ellipticity` to "true" will correct for ellipticity when locating the stations.

- `depth_below_solid_surface`:
If the vertical coordinate of stations is set to "DEPTH", choosing "true" for `depth_below_solid_surface` will use the solid surface as the origin.

- `undulated_geometry`:
If set to "true", the vertical coordinate of stations (either "DEPTH" or "RADIUS") will be determined with respect to the undulated geometry, otherwise with respect to the reference geometry. 

**`wavefields`**

We have now set all the parameters that determine the location of the stations, and we can proceed to set the parameters that will specify which seismic quantities we want to output.

- `coordinate_frame`: Choose from "spz", "RTZ", "ENZ", "xyz". "spz" or (s, phi, z) is the natural reference frame of AxiSEM3D, "RTZ" is (radial, transpose, vertical), "ENZ" is (east, north, vertical), and "xyz" is (x,y,z) in a source-centered frame.

- `medium`: Choose from "SOLID" or "FLUID". `medium` determines whether the stations are contained in a solid or a fluid medium. Note that within a given station group, all stations must belong to the same medium.

- channels:
  
  - Channels available when `medium` set to "SOLID":
    - displacement: U, U1, U2, U3, U_NORM (or |U|)
    - gradient of displacement: G, G11, G12, G13, G21, G22, G23, G31, G32, G33, Gii (or G_I1)
    - strain: E, E11, E12, E13, E21, E22, E23, E31, E32, E33, Eii (or E_I1), E_J2
    - stress: S, S11, S12, S13, S21, S22, S23, S31, S32, S33, Sii (or S_I1), S_J2
    - curl:  R, R1, R2, R3, R_NORM (or |R|)

  - Channels available when `medium` set to "FLUID":
    - displacement: U, U1, U2, U3, U_NORM (or |U|)
    - scalar potential of displacement (U = ∇X): X
    - pressure: P
  - Note that the "1", "2", "3" in the quantities above are determined by the coordinate frame.
  - using U means [U1, U2, U3], and similarly for G, E, S and R; duplicated channels are automatically removed.

**`temporal`**

After having determined where and what we want to output, we need to specify how much of it we want.

- `sampling_period`: Choose from "DT", "DTx2", "DTx3", etc. "DT" stands for the time step in the simulation, hence choosing "DT" will perform no temporal downsampling, "DTx2" will output every other time step, etc. 

- `time_window`: Choose "FULL" to record the whole simulation, or specify a time window [t0,t1].

Finally, all that's left to do is to choose your preferred file options.

**`file_options`**

-`format`: Choose from "ASCII_STATION", "ASCII_CHANNEL", "NETCDF". The first two options outputs the seismograms as a text file. The difference is that "ASCII_STATION" will create individual text files containing all the channels for each station, and is thus limited to a small number of stations; whereas "ASCII_CHANNEL" produces one file per channel, containing all the stations. Finally, "NETCDF" produces a single netcdf file containing everything, and is much more efficient. There is also a parallel version of netcdf which can be activated in CMakeLists.txt.

-`buffer_size`: During the simulation, the solver buffers time steps in RAM, and only writes to file once in a while, because I/O is costly. The drawback is that a larger buffer size increases memory cost, so bear that in mind if you have a lot of stations.

-`flush`: If set to "true", flushes the file after a buffer is written to it. Doing so minimizes data loss in case of abnormal termination of the simulation, but the drawback is that the performance is affected if you chose a small buffer size.

