# inparam.advanced.yaml

It is unlikely that you will need to change anything in this file, as
the options are quite technical and likely of limited interest to the
casual user. 

**geometric deformation**
```yaml
minimum_jacobian: 0.1
```

`minimum_jacobian` is the dimensionless lower limit accepted for the
final particle-relabeling Jacobian after all geometric models have been
applied. The undeformed value is 1 and 0 is singular. If the limit is
violated, the run stops during preloop and reports both the observed and
configured values. Lowering the limit accepts stronger distortion; it
does not smooth or otherwise repair the geometry.

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

`weight_for_load_balancing` may be `ELEMENT_POINT` (using micro-benchmark for METIS in Stage II based on element and point evaluation), `ELEMENT` (same, but ignore point evaluation cost) or `NR` (reuse the same Nr-based weights from Stage-I without timing, this produces a reproducible partition for integration tests).

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
