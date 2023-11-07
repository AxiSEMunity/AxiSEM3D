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
