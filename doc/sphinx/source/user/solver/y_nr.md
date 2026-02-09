# inparam.nr.yaml

This file is only really relevant for when you start including 3D
models. 


For more details, see the section [Nr Field](../nr_field/index.md).

**type**
```
type_Nr: CONSTANT
bound_Nr_by_inplane: true
```

**constant**
```
constant: 5
```

In 1D, in a limited set of cases, and if you are really pushed
for computational resources, you can make a tiny saving by setting
`type_Nr` as `CONSTANT` and equal to `3` (if you’re only using a point
force/dipole source) or `1` if you’re using a monopolar (isotropic
pressure) source instead of the standard `5`.

**analytical**
```
analytical:
    code_ID: depth-dependent (AxiSEM3D default)
    depth_dependent_AxiSEM3D_default:
        control_depths: [0., 50e3, 100e3, 6371e3]
        Nr_at_control_depths: [100, 100, 50, 50]
    any_user_defined_parameters:
        example__bool: true
        example__string: Hello world!
        example__array_of_double: [1., 2., 3.]
        example__array_of_string: [path, file1, file2]
```

**pointwise**
```
pointwise:
    nc_data_file: pointwise.nc
    multip_factor: 1.
```

**structured**
```
structured:
    nc_data_file: structured.nc
    value_out_of_range: 5
```

**wavefield scanning**
```
wavefield_scanning:
    enable_scanning: false
    output_file: scanning_output_Nr.nc
    threshold_Fourier_convergence: 1e-2
    relative_amplitude_skipped: 0.1
    advanced:
        absolute_amplitude_skipped: 1e-12
        max_num_peaks: 10000
        vertex_only: true
        num_steps_per_mesh_period: 12
```



