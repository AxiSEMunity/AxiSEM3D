# inparam.advanced.yaml

## geometric deformation

### `minimum_jacobian`

**What:** minimum acceptable Jacobian after applying geometric models

**Type:** double

**Default:** `0.1`

**Note:**

1) dimensionless; the undeformed value is 1, while 0 is singular
2) checked after all geometric models have been combined; the solver
aborts during preloop if any sampled value is below this threshold
3) must be non-negative; decreasing it permits stronger distortion but
does not repair the input geometry


## verbosity

Parameters for verbosity.

### `channel`

**What:** verbose to stdout or a file

**Type:** string

**Only:** STDOUT / filename

**Default:** `STDOUT`


### `level`

**What:** verbosity level

**Type:** string

**Only:** NONE / ESSENTIAL / DETAILED

**Default:** `ESSENTIAL`


### `warnings`

**What:** show/hide runtime warnings

**Type:** bool

**Default:** `true`


### `loop_info_interval`

**What:** time step interval to display time loop info

**Type:** int

**Default:** `1000`


### `stability_interval`

**What:** time step interval for stability check

**Type:** int

**Default:** `1`

**Note:**

use 1 to precisely locate the origin of instability


## mpi

Parameters for mpi.

### `nproc_per_group`

**What:** number of processors in a mpi group

**Type:** int

**Default:** `1`

**Note:**

1) AxiSEM3D uses a two-level MPI structure where the
processors are divided into groups to avoid broadcasting
a large input dataset (e.g., the Exodus mesh or a large
3D model) on every processor; instead, a large dataset can
be stored only on a "leader" processor in each group, which
handles data requests from its members
2) increase this number (from 1 to the number of processors per
per node) to save memory


### `weight_for_load_balancing`

**What:** weight for load balancing

**Type:** string

**Only:** ELEMENT / ELEMENT_POINT / NR

**Default:** `ELEMENT_POINT`

**Note:**

1) ELEMENT:       use cost measurement on elements
2) ELEMENT_POINT: use cost measurement on both elements and points
3) NR:            notiming, uses Nr(s,z) weights (reproducible)


### `plot_domain_decomposition`

**What:** plot domain decomposition

**Type:** bool

**Default:** `false`

**Note:**

the output netcdf file contains three variables:
* coords,   double, (X, 2), (s,z) of the element centers
* mpi_rank, int,    (X, ),  mpi rank of the elements
* weights,  double, (X, ),  element weights for decomposition
where X is the number of elements


## developers

Parameters for developers.

### `diagnose_preloop`

**What:** enable/disable preloop diagnosis

**Type:** bool

**Default:** `true`

**Note:**

1) output/develop/preloop_diagnosis.log for runtime and memory
2) output/develop/cost_measurements.log for cost measurements
on elements and GLL points


### `max_num_time_steps`

**What:** maximum time steps for running

**Type:** int

**Default:** `0`

**Note:**

use 0 to free this limit


### `time_limit_for_fftw_planning`

**What:** wall-clock time limit (sec) for FFTW planning

**Type:** double

**Default:** `60.`


### `fftw_lucky_numbers`

**What:** enforce FFTW lucky numbers

**Type:** bool

**Default:** `true`

**Note:**

FFTW is good at handling logical sizes of the form:
n = (2^a)*(3^b)*(5^c)*(7^d)*(11^e)*(13^f), where e + f < 2,
as called the lucky numbers; users should use true.


