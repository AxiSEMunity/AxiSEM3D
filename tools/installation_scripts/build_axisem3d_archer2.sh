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
