# inparam.source.yaml







**time axis**

```
time_axis:
    record_length: 1800.
    enforced_dt: NONE
    Courant_number: 1.0
    integrator: NEWMARK
```
 The default
`courant_number` of `1.0` is very unlikely to be stable: start with
something like `0.6` and increase once you know it is working.

<span style="color: red;"><b>*** What is the units of record_length?</b> </span>

----
**sources**

```
list_of_sources:
    - VIRGINIA_201108231751A:
        location:
            latitude_longitude: [37.91, -77.93]
            depth: 12e3
            ellipticity: true
            depth_below_solid_surface: true
            undulated_geometry: true
        mechanism:
            type: MOMENT_TENSOR
            data: [4.71e24, 3.81e22, -4.74e24, 3.99e23, -8.05e23, -1.23e24]
            unit: 1e-7
        source_time_function:
            class_name: GaussianSTF
            half_duration: 50.
            decay_factor: 1.628
            time_shift: 0.
            use_derivative_integral: ERF
            ascii_data_file: stf.txt
            padding: FIRST_LAST
            nc_var_times
            nc_var_data: stf_data
            padding: FIRST_LAST
            chunk_size: NONE
```

The key that you use, e.g., `VIRGINIA_201108231751A` is arbitrary. We normally use the CMT Catalogue designator for the earthquake in question.

We recommend sticking sticking to specfiying the locaiton using `latitude_longitude`, remembering that the north and east is positive and `depth` is in metres – this is important!

The `ellipticity` correction is easy to change and makes little difference unless you need very exact arrival times.

We will come back to `depth_below_solid_surface` and
`undulated_geometry` in the [3D models](../models/index.md) section.

You are most likely to use the `mechanism: type: MOMENT_TENSOR`, in which case you should input the `data` as the components of the moment tensor in the order given by CMTSOLUTION. Remember that CMT uses the non-SI unit of dyne-cm, so you need to add a correction factor by setting `unit` to `1e-7` to turn this into N-m.

For most purposes (seismogram generation) you can do everything involving the `source-time function` (STF) after the simulation is complete in post-processing. In this case, set `half_duration` to `0.0` and the code will make the width of the STF something like 5x the timestep (much less
than the mesh period, i.e. a delta function). Note that in a model like this this is as good as you are going to get – the half-duration cannot be exactly 0.00000, but the difference is miniscule and irrelevant.

Leave `class_name` as `GaussianSTF` and `use_derivative_integral` as either `ERF` or `GAUSSIAN` depending on whether you want displacement or velocity output. We will come back to how to process the output in a later section.

If you are doing something like wavefield visualisation, then you need to set the STF before you do the simulation. In that case, you can set the half-duration and decay factor as you want, or use `netcdf_STF` or
`STREAM_STF` (NetCDF vs. ASCII - more on this later) to read in a particular source-time function that you want to use, for example one downloaded from SCARDEC.

A general note of caution: although you can use any STF you want for seismogram generation, and then deconvolve it out before reconvolving with a new STF, we find that this is rarely a sensible option as it ends up producing numerical artefacts. If you are not sure what STF to use,
use a delta function (as described above) and deal with it later.         