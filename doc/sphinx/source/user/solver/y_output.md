# inparam.output.yaml

**station-wise**
```

list_of_station_groups:
    - global_seismic_network_GSN:
        locations:
            station_file: GSN.txt
            horizontal_x1_x2: LATITUDE_LONGITUDE
            ellipticity: true
            undulated_geometry: true
        wavefields:
            coordinate_frame: RTZ
            medium: SOLID
            channels: [U]
        temporal:
            time_window: FULL
        file_options:
            format: ASCII_STATION
            buffer_size: 1000
            flush: true

```

Here, you need to set the station(s) that you want the displacement or
velocity recorded at. `list_of_station_groups` can take as many networks
(groups of stations) as you want.

We tend to keep `global_seismic_network_GSN` as it is a useful set of
stations all around the world with good data, for testing purposes. If
you are doing something specific though, like array beamforming, you
will want to add your own networks too.

`station_file` is the name of the text file within the input folder
that contains the station coordinates. Subsequent options like
`horizontal_x1_x2` and `vertical_x3` set what the columns in the
`station_file` text file correspond to. For the GSN, these are
`LATITUDE_LONGITUDE` and `DEPTH` (so in your file `GSN.txt`, the columns
are `name, network, latitude, longitude, elevation, depth` – note that
elevation is not used in AxiSEM3D but it is there so that you can swap
files with SPECFEM).

You also need to set the output coordinates using `coordinate_frame` and
*`channels`. We suggest that you leave these as `RTZ` and [*U*] (where
[*U*] = [U1, U2, U3]) for our purposes.

The big choice that you need to make is whether to save your outputs as
text files (`ASCII_STATION` or `ASCII_CHANNEL` depending on what you’re
doing), or netcdf files (netcdf). netcdf is a lot more efficient, but a
little harder to use, especially if you are not used to the file type
before. ASCII has the advantage that you can just open the output and
look at whether it makes sense (if the ground velocity is 10e10 m/s
after 10 seconds, it is probably wrong…) but it is a lot less
space-efficient.

In the processing examples we will stick to ASCII for ease of use but if
you are doing more complex runs or lots of them we would suggest using
NetCDF.

**element-wise**
```
list_of_element_groups: []
```



