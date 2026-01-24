FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    build-essential \
    git \
    curl \
    python3-pip \
    cmake \
    openmpi-bin \
    libopenmpi-dev \
    libeigen3-dev \
    libboost-all-dev \
    libmetis-dev \
    netcdf-bin \
    libnetcdf-dev \
    libfftw3-dev


RUN useradd \
  --create-home \
  axisem3d_user 

USER axisem3d_user

WORKDIR /home/axisem3d_user

RUN mkdir axisem3D_root; cd axisem3D_root; mkdir AxiSEM3D;

COPY . ./axisem3D_root/AxiSEM3D/

RUN cd axisem3D_root; mkdir build; cd build; cmake -Dcxx=mpicxx \
    -Deigen=/usr/include/eigen3 \
    -Dboost=/usr/include/boost \
    -Dfftw=/usr/lib/x86_64-linux-gnu \
    -DFFTW_INCLUDES=/usr/include \
    -Dmetis=/usr/include \
    -Dnetcdf=/usr/include \
    ../AxiSEM3D/SOLVER/

RUN cd axisem3D_root; cd build; make;

ENV PATH="/home/axisem3d_user/axisem3D_root/build:${PATH}"
