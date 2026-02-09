# ASCII file contents

The and define the contents of the ASCII files that we are going to
plot. We assume that you have a file for each station containing three
components of displacement; as set using the `format` and `channels`
options in `inparam.source.yaml`.

Rather than saving the time series (i.e. the list of timesteps) within
each file, this is stored separately in a file called `data_time.ascii`.
This should be a single column of floating point numbers, each of which
is a different point in the time series. The first one is probably
slightly negative, as AxiSEM3D ramps up the source time function before
the zero time index on the simulation. The final entry should be roughly
what you set the record length to be in `inparam.output.yaml`.

If you have not changed the station list, you should have a series of
other files. For example, in the folder `global_seismic_network_GSN/`,
you should have a series of stations, for example `IU.TEIG.ascii`. These
are your station files, and the number of columns in each (and the
order) will depend on the `channels` options you set in
`inparam.output.yaml`. If you used *\[U\]*, you should end up with
three components, in the order R, T, Z if you used `coordinate_frame =
RTZ`.