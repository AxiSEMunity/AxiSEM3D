# List of Element Groups

## elements to be recorded

### `horizontal_range`

**What:** horizontal range

**Type:** array of double

**Default:** `[0, 1e10]`

**Note:**

1) θ-range for a spherical mesh
2) s-range for a Cartesian mesh


### `vertical_range`

**What:** vertical range

**Type:** array of double

**Default:** `[0, 1e10]`

**Note:**

1) r-range for a spherical mesh
2) z-range for a Cartesian mesh


## in-plane sampling on each element

Parameters for in-plane sampling on each element.

### `edge_dimension`

**What:** edge dimension

**Type:** string

**Only:** HORIZONTAL / VERTICAL / BOTH

**Default:** `BOTH`

**Note:**

1) use HORIZONTAL or VERTICAL to record only an outer edge
of each element, with the edge position specified below
2) use BOTH to record the entire element


### `edge_position`

**What:** edge position

**Type:** double

**Default:** `6371e3`

**Note:**

1) as we only search among the four outer edges of each
element, this edge position must lie on a mesh
boundary (either internal or external), such as one
of the vertical discontinuities
2) not used if edge_dimension = BOTH


### `GLL_points_one_edge`

**What:** GLL points on one edge

**Type:** string / array of int

**Only:** FULL / a subset of [0, 1, 2, ..., npol], where npol is

**Default:** `[0,4]`

**Note:**

use this parameter for in-plane downsampling;
for example, when npol = 4:
- [2]:       only the element center
- [0, 4]:    only the 4 vertices
- [0, 2, 4]: 4 vertices + 4 edge centers + 1 element center
- FULL:      all GLL points, the same as [0, 1, 2, 3, 4]


## azimuthal sampling on each element

### `phi_list`

**What:** azimuths (ɸ) of the recorded slices

**Type:** array of double

**Default:** `[0, 0.838]`

**Note:**

the unit is radian
phi_list: [0, 0.838, 1.57079632679, 2.409, 3.14159265359, 3.979, 4.71238898038, 5.55]


### `lat_lon_list`

**What:** latitudes and longitudes of the recorded slices

**Type:** array of array of double

**Default:** `[]`

**Note:**

ɸ computed from (lat, lon) will be append to phi_list


### `na_space`

**What:** grid space for azimuthal storage

**Type:** int

**Default:** `1`

**Note:**

1) when the Fourier series coefficients are recorded
(phi_list = [] and lat_lon_list = []), this parameter
affects how the wavefield data with different azimuthal
dimensions are stored in the output NetCDF file(s);
see the online manual for details
2) must be greater than or equal to 1
3) recommended: maximum Nr in the target region / 200


## wavefields to be recorded

Parameters for the wavefields to be recorded.

### `coordinate_frame`

**What:** coordinate frame of wavefields

**Type:** string

**Only:** spz / RTZ / ENZ

**Default:** `RTZ`

**Note:**

1) spz: (s, phi, z) or AxiSEM3D-intrinsic
2) RTZ: (radial, transpose, vertical)
3) ENZ: (east, north, vertical)


### `medium`

**What:** type of medium

**Type:** string

**Only:** SOLID / FLUID

**Default:** `SOLID`

**Note:**

all elements in a group must be located in either
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

**Default:** `0.1`

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

### `buffer_size`

**What:** number of sampled time steps to be buffered

**Type:** int

**Default:** `100`

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


## elements to be recorded

### `horizontal_range`

**What:** horizontal range

**Type:** array of double

**Default:** `[0, 1e10]`

**Note:**

1) θ-range for a spherical mesh
2) s-range for a Cartesian mesh


### `vertical_range`

**What:** vertical range

**Type:** array of double

**Default:** `[0, 1e10]`

**Note:**

1) r-range for a spherical mesh
2) z-range for a Cartesian mesh


## inplane sampling on each element

Parameters for inplane sampling.

### `edge_dimension`

**What:** edge dimension

**Type:** string

**Only:** HORIZONTAL / VERTICAL / BOTH

**Default:** `VERTICAL`

**Note:**

1) use HORIZONTAL or VERTICAL to record only an outer edge
of each element, with the edge position specified below
2) use BOTH to record the entire element


### `edge_position`

**What:** edge position

**Type:** double

**Default:** `6371.1e3`

**Note:**

1) as we only search among the four outer edges of each
element, this edge position must lie on a mesh
boundary (either internal or external), such as one
of the vertical discontinuities
2) not used if edge_dimension = BOTH


### `GLL_points_one_edge`

**What:** GLL points on one edge

**Type:** string / array of int

**Only:** FULL / a subset of [0, 1, 2, ..., npol], where npol is

**Default:** `[0, 2, 4]`

**Note:**

use this parameter for inplane downsampling;
for example, when npol = 4:
- [2]:       only the element center
- [0, 4]:    only the 4 vertices
- [0, 2, 4]: 4 vertices + 4 edge centers + 1 element center
- FULL:      all GLL points, the same as [0, 1, 2, 3, 4]


## azimuthal sampling on each element

### `phi_list`

**What:** azimuths (ɸ) of the recorded slices

**Type:** array of double

**Default:** `[]`

**Note:**

the unit is radian


### `lat_lon_list`

**What:** latitudes and longitudes of the recorded slices

**Type:** array of array of double

**Default:** `[]`

**Note:**

ɸ computed from (lat, lon) will be append to phi_list


### `na_space`

**What:** grid space for azimuthal storage

**Type:** int

**Default:** `1`

**Note:**

1) when the Fourier series coefficients are recorded
(phi_list = [] and lat_lon_list = []), this parameter
affects how the wavefield data with different azimuthal
dimensions are stored in the output NetCDF file(s);
see the online manual for details
2) must be greater than or equal to 1
3) recommended: maximum Nr in the target region / 200


## wavefields to be recorded

Parameters for wavefields to be recorded.

### `coordinate_frame`

**What:** coordinate frame of wavefields

**Type:** string

**Only:** spz / RTZ / ENZ

**Default:** `RTZ`

**Note:**

1) spz: (s, phi, z) or AxiSEM3D-intrinsic
2) RTZ: (radial, transpose, vertical)
3) ENZ: (east, north, vertical)


### `medium`

**What:** type of medium

**Type:** string

**Only:** SOLID / FLUID

**Default:** `SOLID`

**Note:**

all elements in a group must be located in either
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

**Default:** `0.1`

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

### `buffer_size`

**What:** number of sampled time steps to be buffered

**Type:** int

**Default:** `100`

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


