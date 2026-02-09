
# The Code
The mesher is written in Python, and the solver in modern C++ while relying on MPI parallelisation. The solver requires netcdf, fftw, boost and eigen libraries. It interfaces with common seismological I/O formats and libraries such as CMTSOLUTION, SAC, ObsPy, ASDF. The code runs very efficiently on parallel infrastructures and has shown excellent scalability up to 12,000 cores and more. 


 ## Languages 

The AxiSEM3D code (Lent et al., 2016;2019) is written in the **C++** programming language.
Various sub-routines involve Fortran, and some of the pre- and
post-processing codes are supplied in Python or Matlab. Input files are
human-readable text files, and no particular programming languages are
required to understand and use them.

The casual user will not need to edit the source code, and as such, no
familiarity with C++ is required. You may see references elsewhere to
the ‘old’ AxiSEM code (Nissen-Meuyer et al., 2014), which is written in Fortran. You do not need
to be familiar with either old AxiSEM or Fortran to run AxiSEM3D.

## Understanding inputs and outputs 

The file formats associated with inputs and outputs are (at least
compared to the rest of the code) simple. Unless you choose otherwise,
you should be able to read these as plain text.

For viewing and editing input files, you may wish to consider using an
editor such as **VSCode** or **Kate** *according to Jonathan, Kate is
amazing; no one else has ever heard of her*, which will render the text
in the files in an easier-to-read way (for example, keys and values in
different colours). Else, readers such as vim or gedit are also fine.

At the simplest level, you can save outputs as ASCII files and plot them
using something as simple as Excel. As discussed in the user guide,
non-human readable formats are more space efficient, but do require some
basic familiarity with Python or Matlab (or equivalent) to read and
plot.

If you are not biased either way, we suggest doing post-processing in
**Python**, since you can then make use of the enormous functionality
offered by the **ObsPy** (Beyreuther et al., 2010) package.

## AxiSEM3D architecture 

Upon downloading the source code from Github (more on this later), it is
worth taking a brief look at it, even if you do not plan on doing any
coding in C++ yourself.

You will see that the main body of the code is divided into two
sections: the **core** and the **preloop**.

### The Preloop 

The preloop involves everything that only needs to be done once in a
simulation before it starts – e.g., dividing up the mesh between
different nodes. The code will also perform a series of checks in the
preloop, for example checking that the Jacobian for each element (used
to translate between the physical configuration and the reference
configuration used for numerical integration) is positive.

### The Core 

The core of the code contains the routines for executing the main time
loop, that is, solving the equations of motion numerically at each
timestep.

AxiSEM3D is designed to be highly parallel, meaning that the workload in
the time loop is efficiently divided between a number of different
computational nodes and their constituent processors. The communication
overheads between different nodes, and the time required to assemble the
whole solution at the end of the simulation, are designed to be small
compared to the overall runtime of the code.

## Computing architectures 

The vocabulary surrounding computing architectures can be difficult to
get to grips with at first. Here, we will briefly discuss the
relationship between processors and nodes, as applied to AxiSEM3D.

### On your machine

When running an AxiSEM3D simulation on a local machine (e.g., your
laptop), you generally only need to set the number of **processes** that
you want to run. In short, a process is a parallel ‘thread’ of
operations that is executed simultaneously (‘in parallel’) to other
processes. Most of the time, this should be the same as the number of
processors that your machine has, though some machines support more
advanced operations such as hyperthreading. Generally, creating more
processes than you have processors starts to reduce efficiency again.

### On an HPC architecture 

On **high-performance computing** (HPC) architectures, i.e. **clusters
and supercomputers**, you will probably need to be a bit more specific
about how you want the code to execute. In general, you will need to
specify the number of nodes, the number of processors, and the runtime.

Most HPC systems have a standard number of processors per node, say 64
or 128. This means that you can have up to this many independent threads
running on a single shared-memory unit. Each node will also have a fixed
amount of memory (around 256 or 512 GB often) which all the processors
will have to share between them.

Sometimes, you will also have the option to request **high-memory
nodes**, which can be more efficient if you are doing memory-intensive
computations (specifically, things like wavefield visualisation which
have a heavy input/output load). These high-memory nodes can let you run
a higher number of processors-per-node without getting an Out Of Memory
(OOM) error than the equivalent standard nodes; unused processors on a
node are left dormant and hence the simulation is less efficient if this
occurs.

Note that most HPC systems have separate **login and compute nodes**. It
is important to make sure that you only run AxiSEM3D on the compute
nodes, as these have the required power and memory to execute
simulations. Running on the login nodes will likely result in error from
a lack of permissions, slowing down the machine for all users, or even
crashing it.