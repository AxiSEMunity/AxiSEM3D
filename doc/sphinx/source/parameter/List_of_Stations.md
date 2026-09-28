# List of Stations

## station location

Parameters for station location.

### `station_file`

**What:** station location file

**Type:** filename

**Default:** `STATIONS_OUTPUT.txt`

**Note:**

1) must be an ascii file with five or six columns:
NAME NETWORK x1 x2 [useless] x3
[useless] is for compatibility with SPECFEM and AxiSEM
2) (x1, x2, x3) are specified by the next two options
3) empty lines and comment lines (led by #) are allowed


### `horizontal_x1_x2`

**What:** horizontal coordinates x1 and x2

**Type:** string

**Only:** LATITUDE_LONGITUDE / DISTANCE_AZIMUTH / XY_CARTESIAN

**Default:** `XY_CARTESIAN`

**Note:**

1) the unit for LATITUDE_LONGITUDE is degree
2) the unit for DISTANCE is either radian (for a
spherical mesh) or meter (for a Cartesian mesh);
the unit for AZIMUTH is radian
3) the unit for XY_CARTESIAN is meter


### `vertical_x3`

**What:** vertical coordinate x3

**Type:** string

**Only:** RADIUS / DEPTH

**Default:** `DEPTH`

**Note:**

the unit is meter


### `ellipticity`

**What:** correct for ellipticity when locating the stations

**Type:** bool

**Default:** `false`

**Note:**

used only when horizontal_x1_x2 = LATITUDE_LONGITUDE


### `depth_below_solid_surface`

**What:** use solid surface as depth origin

**Type:** bool

**Default:** `true`

**Note:**

used only when vertical_x3 = DEPTH


### `undulated_geometry`

**What:** use undulated (otherwise reference) geometry to

**Type:** bool

**Default:** `true`

**Note:**

compatible with vertical_x3 = both RADIUS and DEPTH


## wavefields to be recorded

Parameters for wavefields.

### `coordinate_frame`

**What:** coordinate frame of wavefields

**Type:** string

**Only:** spz / RTZ / ENZ / xyz

**Default:** `ENZ`

**Note:**

1) spz: (s, phi, z) or AxiSEM3D-intrinsic
2) RTZ: (radial, transpose, vertical)
3) ENZ: (east, north, vertical)
4) xyz: (x, y, z) in source-centered frame


### `medium`

**What:** type of medium

**Type:** string

**Only:** SOLID / FLUID

**Default:** `SOLID`

**Note:**

all stations in a group must be located in either
the solid or the fluid domain


### `channels`

**What:** list of channels

**Type:** array of string

**Default:** `[U]`

**Note:**

1) allowed channels for medium = SOLID
* displacement:
U, U1, U2, U3, U_NORM (or |U|)
* gradient of displacement:
G, G11, G12, G13, G21, G22, G23, G31, G32, G33,
Gii (or G_I1)
* strain:
E, E11, E12, E13, E21, E22, E23, E31, E32, E33,
Eii (or E_I1), E_J2
* stress:
S, S11, S12, S13, S21, S22, S23, S31, S32, S33,
Sii (or S_I1), S_J2
* curl:
R, R1, R2, R3, R_NORM (or |R|)
2) allowed channels for medium = FLUID
* displacement:
U, U1, U2, U3, U_NORM (or |U|)
* scalar potential of displacement (U = ∇X):
X
* pressure:
P
3) (1, 2, 3) are determined by coordinate_frame
4) using U means [U1, U2, U3], and similarly for G, E, S
and R; duplicated channels are automatically removed


## temporal sampling

Parameters for temporal sampling.

### `sampling_period`

**What:** sampling period

**Type:** string / double

**Only:** DT / DTx2 / DTx3 / ... / value

**Default:** `0.2`

**Note:**

DT stands for Δt of the simulation; DTx3 means
sampling period = Δt * 3


### `time_window`

**What:** time window

**Type:** string / array of double

**Only:** FULL / [t0, t1]

**Default:** `FULL`

**Note:**

use FULL to record the whole simulation


## file options

Parameters for file options.

### `format`

**What:** output file format

**Type:** string

**Only:** ASCII_STATION / ASCII_CHANNEL / NETCDF

**Default:** `NETCDF`

**Note:**

1) ASCII_STATION: one ascii file contains all channels at
one station, available only for a small
number of stations
2) ASCII_CHANNEL: one ascii file contains one channel at
all stations, available for many stations
3) NETCDF: much more efficient than ascii, available for
many stations; parallel NetCDF can be activated
in CMakeLists.txt


### `buffer_size`

**What:** number of sampled time steps to be buffered

**Type:** int

**Default:** `1000`

**Note:**

1) the solver buffers wave data during the time loop for
efficient writing; increase this buffer size to save
writing time and decrease it to save memory (useful for
a large number of stations)
2) this parameter does not affect the final results


### `flush`

**What:** flush file after writing a buffer to it

**Type:** bool

**Default:** `true`

**Note:**

1) pro: minimizes data loss in case of abnormal termination
con: hits output performance if buffer_size is small
2) this parameter does not affect the final results


## element-wise

### `list_of_element_groups`

**What:** list of element groups

**Type:** array of objects

**Default:** `[]`

**Note:**

1) different options such as channels and sampling rates can
be used for different element groups
2) use [] if no element group presents


