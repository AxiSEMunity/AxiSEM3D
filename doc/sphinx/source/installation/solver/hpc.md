# HPC Cluster Installation

Most HPC clusters already provide these libraries through environment modules.

## Browse Available Modules

```bash
module avail
module avail fftw
```

## Load Required Packages

```bash
module load fftw
module load fftw/3.3.8
```

If a required library is missing, you may:

1. contact your system administrator  
2. install dependencies locally via conda (if permitted)

> **Note:** Modules loaded during compilation are usually required at runtime as well.

---

## Configure the Build on HPC

On most clusters, dependencies are provided via environment modules.  
After loading the required modules, installation prefixes are usually available
as environment variables such as `HDF5_ROOT`, `NETCDF_ROOT`, etc.

For example:

```bash
module load hdf5 netcdf fftw metis boost eigen
module show hdf5
echo $HDF5_ROOT
```

You can then configure the build using these variables directly:

```bash
rm -rf build && cmake -S ./SOLVER -B build \
  -Dcxx=mpicxx \
  -Dhdf5=$HDF5_ROOT \
  -Dnetcdf=$NETCDF_ROOT \
  -Deigen=$EIGEN3_ROOT \
  -Dboost=$BOOST_ROOT \
  -Dfftw=$FFTW_ROOT \
  -Dmetis=$METIS_ROOT
```

> If your cluster does not provide these variables, you can inspect module paths via `module show <package>`.

---

## Compile

```bash
cmake --build build -j4
```

### Testing in parallel

Most HPC systems do not allow you to just use resources at will. We
assume that you have an account on the system that you are using, and
that there are some resources attached to that account.

In order to run in parallel, you need to request a certain number of
nodes, and a runtime. On most machines, you also need to tell the system
how many processors you want to run on each node. Remember, the total
memory of each node is fixed and must be partitioned between all
processors, so for testing purposes it is sensible to set the number of
computing processors-per-node to the actual number of
processors-per-node, as any leftover processors cannot be assigned to
other users and are just wasted.

### Interactive mode 

The best way to test and debug the code on HPC systems is in
**interactive mode**. Interactive mode is different to **batch/schedule
mode**. In the former, you tell the system that you need X nodes for Y
minutes (a sensible test case might be 4 nodes, each with 128 processors
for example, for 10 minutes at a time). The system will generally
allocate them to you immediately, or after a very short delay.

In the latter, you do the same, but the request can normally be larger
(sometimes a lot larger - more on this later). The delay for this many
processors to become available can be significant, often hours to days
on busy machines for large jobs. So, testing in interactive mode is much
faster.

Once you have been assigned your interactive mode session, you will
probably need to re-activate your Conda environment, etc. You can do
this in the same way that you did it before. You do not need to re-build
the code, but if you loaded any modules ‘by hand’, i.e. using something
like module load fftw in the command line you will need to do this
again. To save time, we normally use a configuration script to do this
(more on this later).

Then, to test in interactive mode, you will need to use **mpirun**, or
sometimes an alternative called **aprun**, to run in parallel.

You can do this by typing **mpirun -np X ./build/axisem3d**, or **aprun -n X
./build/axisem3d**. Note that X is the number of processes to spawn (in this
case, equal to the number of processors that the code will run on), not
the number of nodes. The most efficient solution is to set X to the
number of processors-per-node multiplied by the number of nodes that you
requested in the interactive session.
