# Advanced settings

The input file `inparam.advanced.yaml` contains several advanced settings. Most of these parameters facilitate code development and are usually not of users' interest except `nproc_per_group` under `mpi`, which is important for large-scale computations.


**`verbose`**

The parameters under `verbose` control the verbosity of the standard output channel. The only parameter the user is likely to change is `channel`, which can direct the output to a file; an equivalent way of doing this is 

```bash
mpirun -np 4 axisem3d > std_output.txt
```

**`mpi`**

The parameters under `mpi` controls the MPI behavior of the solver, among which `nproc_per_group` is the most useful one. Some data are required by all processors, such as the material parameters of a 3D model. We can broadcast and store such data on every processor; however, because these data are not distributed, the memory usage will be too large. Alternatively, we can store them only on the `ROOT` processor, which handle the request from the other processors. This way, we minimize the memory, but the runtime (e.g., for reading a large 3D model) will be too long, as the `ROOT` has to process the requests from the others sequentially. To solve this conflict, we adopt a two-level MPI architecture. The data will be broadcasted and stored on a small group of processors, which we call the *leader* processors, which handle requests of its *followers*. The parameter `nproc_per_group` determines how many followers are managed by one leader. 

How to choose `nproc_per_group`:

* the smallest value is 1, which is most time-efficient and memory-costly;
* the largest value is the number of processors per node, which varies with the HPC you use, which is most memory-efficient and slow.

Any number in between is reasonable. If the program is terminated by a memory error, increase this number. If the memory error persists if you have reached the number of processors per node, it means the number of nodes is not enough for the problem you are trying to solve.

**`develop`**
The parameters under `develop` are mostly for developers. `max_num_time_steps` can be a useful one for users, which allows you to run only one time step to validate the input files. 

