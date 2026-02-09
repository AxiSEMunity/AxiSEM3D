# inparam.advanced.yaml

It is unlikely that you will need to change anything in this file, as
the options are quite technical and likely of limited interest to the
casual user. 

**verbosity**
```
verbose:
    channel: STDOUT
    level: ESSENTIAL
    warnings: true
    loop_info_interval: 1000
    stability_interval: 1
```
**mpi**
```
mpi:
    nproc_per_group: 1
    weight_for_load_balancing: ELEMENT_POINT
    plot_domain_decomposition: false
```

**developers**
```
develop:
    diagnose_preloop: true
    max_num_time_steps: 0
    time_limit_for_fftw_planning: 60.
    fftw_lucky_numbers: true
```
Update `max_num_time_steps` if you want to make sure
that the code will at least run through a small number of timesteps, e.g., 10,
without becoming unstable.