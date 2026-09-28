# List of 3D Models

## modelname

Parameters for the model.

### `activated`

**What:** activate this model

**Type:** bool

**Default:** `true`


### `class_name`

**What:** class name

**Type:** string

**Default:** `StructuredGridV3D`

**Note:**

current built-in classes include
- StructuredGridV3D: volumetric 3D model on a structured grid
- StructuredGridG3D: geometric 3D model on a structured grid
- StructuredGridO3D: ocean-load 3D model on a structured grid
- Ellipticity: deform the mesh with global ellipticity


### `nc_data_file`

**What:** NetCDF data file

**Type:** filename

**Default:** `lowermost_mantle_ani.nc`


### `horizontal`

**What:** type of horizontal coordinates

**Type:** string

**Only:** DISTANCE_AZIMUTH / XY_CARTESIAN / LATITUDE_LONGITUDE

**Default:** `LATITUDE_LONGITUDE`


### `vertical`

**What:** type of vertical coordinate

**Type:** string

**Only:** RADIUS / DEPTH

**Default:** `DEPTH`


### `ellipticity`

**What:** correct for ellipticity when locating the model

**Type:** bool

**Default:** `false`

**Note:**

used only when horizontal = LATITUDE_LONGITUDE


### `depth_below_solid_surface`

**What:** use solid surface as depth origin

**Type:** bool

**Default:** `false`

**Note:**

used only when vertical = DEPTH


### `nc_variables`

**What:** NetCDF variables for the coordinates

**Type:** array of string

**Default:** `[latitude, longitude, depth]`


### `data_rank`

**What:** rank of the coordinates in data

**Type:** array of int

**Default:** `[1, 2, 0]`


### `length_unit`

**What:** length unit of the coordinates

**Type:** string / value

**Only:** km / m / number

**Default:** `km`


### `angle_unit`

**What:** angle unit of the coordinates

**Type:** string

**Only:** degree / radian

**Default:** `degree`


### `undulated_geometry`

**What:** use undulated (otherwise reference) geometry to

**Type:** bool

**Default:** `false`

**Note:**

compatible only with vertical = RADIUS


### `whole_element_inplane`

**What:** check inplane model range for the whole element

**Type:** bool

**Default:** `false`

**Note:**

1) if this parameter is set to true, the element center
will be used to determine whether an element is located
within the "inplane" model range
2) if its center is in range, all its GLL points must be
in range, or an exception will occur; users can extend
the model range slightly to allow for numerical errors
3) this parameter safely realizes inplane discontinuities


## parameters for properties

Parameters for model properties.

### `nc_var`

**What:** NetCDF variable

**Default:** `C11`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C12`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C13`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C14`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C15`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C16`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C22`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C23`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C24`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C25`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C26`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C33`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C34`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C35`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C36`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C44`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C45`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C46`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C55`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C56`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `nc_var`

**What:** NetCDF variable

**Default:** `C66`


### `factor`

**What:** factor or unit

**Default:** `1e9`


### `reference_kind`

**What:** reference kind

**Only:** ABS / REF1D / REF3D / REF_PERTURB

**Default:** `ABS`

**Note:**

For any property X:
1) ABS: absolute value
X_3D = value_in_file
2) REF1D: perturbation w.r.t. the 1D reference model
X_3D = (1 + value_in_file) * X_1D
3) REF3D: perturbation w.r.t. the current 3D model
X_3D = (1 + value_in_file) * X_3D
4) REF_PERTURB => perturbation w.r.t. the current
perturbation or (X_3D - X_1D)
X_3D = (1 + value_in_file) * (X_3D - X_1D) + X_1D


### `store_grid_only_on_leaders`

**What:** store grid data only on the leader processors

**Type:** bool

**Default:** `true`

**Note:**

turn this on if the model is large; set mpi:nproc_per_group
in inparam.advanced.yaml to the number of processors per
node to minimize memory usage


