# List of Sources

## source location

Parameters for source location.

### `latitude_longitude`

**What:** horizontal location

**Type:** string / array of double

**Only:** ON_AXIS / values

**Default:** `[89.99999,0.0]`

**Note:**

1) the key can be "latitude_longitude" or "distance_azimuth"
2) the unit for "latitude_longitude" is degree;
the unit for "distance" is either radian (for a
spherical mesh) or meter (for a Cartesian mesh);
the unit for "azimuth" is radian
3) use ON_AXIS to place the source on the mesh axis
4) if geodesy::lat_lon_north_pole_mesh = SOURCE in
inparam.model.yaml and this is the FIRST source in this
list, the key must be "latitude_longitude"


### `depth`

**What:** vertical location

**Type:** double

**Default:** `0.0`

**Note:**

the key can be "depth" or "radius" and the unit is meter


### `ellipticity`

**What:** correct for ellipticity when locating the source

**Type:** bool

**Default:** `true`

**Note:**

used only when the key for horizontal location is
"latitude_longitude"


### `depth_below_solid_surface`

**What:** use solid surface as depth origin

**Type:** bool

**Default:** `true`

**Note:**

used only when the key for vertical location is "depth"


### `undulated_geometry`

**What:** use undulated (otherwise reference) geometry to

**Type:** bool

**Default:** `true`

**Note:**

compatible with both "depth" and "radius"


## source mechanism

Parameters for source mechanism.

### `type`

**What:** type of source mechanism

**Type:** string

**Only:** MOMENT_TENSOR, FORCE_VECTOR, FLUID_PRESSURE

**Default:** `MOMENT_TENSOR`


### `data`

**What:** data for the source mechanism

**Type:** array of double

**Default:** `[2e25, -1e25, -1e25, 0, 0, 0]`

**Note:**

1) use [M11, M22, M33, M12, M13, M23] for MOMENT_TENSOR;
[F1, F2, F3] for FORCE_VECTOR;
[P] for FLUID_PRESSURE,
where 123 stands for ZRT (vertical, radial, transpose)
2) if horizontal location is given by "latitude_longitude",
the RT-axes are determined w.r.t. the north pole;
the moment tensor of an earthquake then follows the same
order as in the CMTSOLUTION format (globalcmt.org)
3) if horizontal location is given by "distance_azimuth",
the RT-axes are determined w.r.t. the source (mesh axis)


### `unit`

**What:** unit of data

**Type:** double

**Default:** `1`

**Note:**

use 1e-7 to convert dyn*cm (in CMTSOLUTION) to N*m


## source time function

Parameters for source time function.

### `class_name`

**What:** class name

**Type:** string

**Only:** GaussianSTF / StreamSTF / NetCDF_STF

**Default:** `GaussianSTF`


### `half_duration`

**What:** half duration

**Type:** double

**Default:** `0.0`

**Note:**

use 0 for a Dirac-δ to compute the Green's functions


### `decay_factor`

**What:** decay factor

**Type:** double

**Default:** `1.628`

**Note:**

1) this factor further controls the shape of the bell
2) SPECFEM3D_GLOBE uses 1.628 and AxiSEM 3.5


### `time_shift`

**What:** time shift

**Type:** double

**Default:** `0.`

**Note:**

1) the axis of the bell will be aligned at the time
origin (t=0) if time_shift is given by zero
2) a positive time shift moves the bell to the right


### `use_derivative_integral`

**What:** use time derivative or integral of Gaussian

**Type:** string

**Only:** ERF / GAUSSIAN / FIRST_DERIVATIVE / RICKER

**Default:** `GAUSSIAN`

**Note:**

1) ERF: error function, integral of Gaussian
2) GAUSSIAN: Gaussian
3) FIRST_DERIVATIVE: first derivative of Gaussian
4) RICKER: second derivative of Gaussian


### `ascii_data_file`

**What:** ascii data file

**Type:** filename

**Default:** `stf.txt`

**Note:**

1) this file must contain two columns, the first
for times and the second for values, and at least
two rows (two time-value points)
3) empty lines and comment lines (led by #) are allowed


### `padding`

**What:** padding values when time is out of range

**Type:** string / array of double

**Only:** NONE / FIRST_LAST / [left, right]

**Default:** `FIRST_LAST`

**Note:**

1) NONE: throw an error when time is out of range
2) FIRST_LAST: use the first and the last values in file
3) [left, right]: specified padding values


### `nc_data_file`

**What:** netcdf data file

**Type:** filename

**Default:** `stf_file.nc`


### `nc_var_times`

**What:** netcdf variable for times

**Type:** string

**Default:** `stf_times`

**Note:**

must be double precision with shape (X,), X>=2


### `nc_var_data`

**What:** netcdf variable for data

**Type:** string

**Default:** `stf_data`

**Note:**

1) must have the same shape as nc_var_times
2) datatype must match solver precision
(refer to USE_DOUBLE in CMakeLists.txt)


### `padding`

**What:** padding values when time is out of range

**Type:** string / array of double

**Only:** NONE / FIRST_LAST / [left, right]

**Default:** `FIRST_LAST`

**Note:**

1) NONE: throw an error when time is out of range
2) FIRST_LAST: use the first and the last values in file
3) [left, right]: specified padding values


### `chunk_size`

**What:** chunk size (how many time steps per chunk)

**Type:** NONE / int

**Default:** `NONE`

**Note:**

1) NetCDF_STF loads data by chunks during the time loop;
increase this chunk size to save reading time and
decrease it to save memory (useful for a large number
of sources); must be greater than 2.
2) use NONE to disable buffering (chunk size = file size)


