# Example Installations on HPC Systems
This page provides details of our experience compiling AxiSEM-3D on a number of high
performance computing facilities. This is often an interesting experience and, because
the systems and AxiSEM-3D are both evolving there is no guarantee that the solution we
found still works now. Nevertheless, the information below may be of some help to
others.

## Archer 2 

We can mostly compile AxiSEM-3D using libraries provided on the system. However,
the latest boost is (one minor version) too old so that needs to be downloaded.
There is also an interesting issue with the naming scheme of metis, which is
avoided by making a local soft link to the library with a more common name
(otherwise CMake complains). Finally, it is important to use the gcc compiler
as the Cray C++ compiler throws an ICE which we have not investigated.

These commands worked on 8-4-2022:

```
#!/bin/bash

PRFX=/home/n03/n03/andreww/AxiSEM-3D

cd ${PRFX}

module load PrgEnv-gnu
module load cray-fftw/3.3.8.9
module load metis/5.1.0
module load cray-hdf5/1.12.0.3
module load cray-netcdf/4.7.4.3
module load cmake/3.21.3

mkdir -p dependencies
cd dependencies

[ ! -d ./boost_1_78_0 ] && wget -c https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.bz2 -O - | tar -jx -C ./

mkdir -p ./metis/lib
ln -s /work/y07/shared/libs/core/metis/5.1.0/GNU/9.3/bin ./metis/bin
ln -s /work/y07/shared/libs/core/metis/5.1.0/GNU/9.3/include ./metis/include
ln -s /work/y07/shared/libs/core/metis/5.1.0/GNU/9.3/lib/libmetis_gnu.a ./metis/lib/libmetis.a

cd ${PRFX}

git clone https://github.com/kuangdai/AxiSEM-3D.git AxiSEM3D
git -C AxiSEM3D pull

mkdir -p build
cd build

export CRAYPE_LINK_TYPE=dynamic

rm -rf ./*

cmake -Dcxx=CC \
  -Dflags="-std=c++17 -O3 -DNDEBUG -fPIC" \
  -Deigen=/work/y07/shared/libs/core/eigen/3.4.0/include \
  -Dboost=${PRFX}/dependencies/boost_1_78_0 \
  -Dfftw=/opt/cray/pe/fftw/3.3.8.9/x86_rome \
  -Dmetis=${PRFX}/dependencies/metis \
  -Dnetcdf=/opt/cray/pe/netcdf/4.7.4.3/GNU/9.1 \
  -Dhdf5=/opt/cray/pe/hdf5/1.12.0.3/GNU/9.1 \
  -Dpar_netcdf=false \
  -Dnpol=4 \
  -Ddouble=true \
  ../AxiSEM3D/SOLVER/

make -j
```
  
## ARC (University of Oxford)

The 'trick' here is to mix up-to-date header only C++ libraries and link to 
compiled libraries from the system (using the module system and letting
CMake find the link paths without our 'help'). The default compiler is
g++ (this is what mpicxx points to) and that seems to be the sensible
choice. Anyway, as of  9/4/2022 this works:

```
#!/bin/bash

PRFX=/home/eart0526/code/AxiSEM3D_april

cd ${PRFX}

module load netCDF
module load CMake

# mpicxx points to g++, which should be okay I think

mkdir -p dependencies
cd dependencies

[ ! -d ./eigen-3.4.0 ] && \
        wget -c https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2 -O - | tar -jx -C ./

[ ! -d ./boost_1_78_0 ] && \
	wget -c https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.bz2 -O - | tar -jx -C ./


cd ${PRFX}

git clone https://github.com/kuangdai/AxiSEM-3D.git AxiSEM3D
git -C AxiSEM3D pull

mkdir -p build
cd build


rm -rf ./*

cmake -Dcxx=mpicxx \
  -Dflags="-std=c++17 -O3 -DNDEBUG -fPIC" \
  -Deigen=${PRFX}/dependencies/eigen-3.4.0 \
  -Dboost=${PRFX}/dependencies/boost_1_78_0 \
  -Dfftw=/apps/system/easybuild/software/FFTW/3.3.10-gompi-2021b \
  -Dmetis=/apps/system/easybuild/software/METIS/5.1.0-GCCcore-10.3.0 \
  -Dpar_netcdf=false \
  -Dnpol=4 \
  -Ddouble=true \
  ../AxiSEM3D/SOLVER/

make
```

There are a couple of system specific things to note. First, the login nodes
are low-powered VMs with a different architecture to the compute nodes. This
means sensible building must be done on an interactive node requested via slurm
with `srun -p interactive --pty /bin/bash`. That just requests one core, but
`make -j` will try and run on all cores and end up getting killed mid-build.
This is why the `-j` flag is omitted above. Just get a coffee.

Finally, you'll need to load modules at run time (in the slurm submission
script). This works for me:

```
#!/bin/bash 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=48
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:10:00 
#SBATCH --job-name=myjob 
#SBATCH --partition=devel


module load netCDF
module load FFTW

# Why does it expect input based on the path to the exe??
cp /home/eart0526/code/AxiSEM3D_april/build/axisem3d .
mpirun axisem3d 
rm axisem3d
```

## GRACE cluster at Yale

AxiSEM3D is installed on GRACE as a module. Simply load it via module
load AxiSEM3D/2022Sep13-iomkl-2020b-avx512. To run it, just use mpirun
axisem3d.
