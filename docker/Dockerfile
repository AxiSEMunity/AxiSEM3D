FROM ubuntu:questing

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
    libfftw3-dev \
    libhdf5-dev


RUN useradd \
  --create-home \
  axisem

USER axisem

WORKDIR /home/axisem/

RUN git clone https://github.com/AxiSEMunity/AxiSEM3D.git AxiSEM3D;

RUN cd AxiSEM3D; mkdir build; cd build; cmake -B . \
    -Deigen=/usr/include/eigen3 \
    -Dboost=/usr/include/boost \
    -Dfftw=/usr/lib/x86_64-linux-gnu \
    -DFFTW_INCLUDES=/usr/include \
    -Dmetis=/usr/include \
    -Dnetcdf=/usr/include \
    -Dhdf5=/usr \
    ../

#/usr/lib/x86_64-linux-gnu/hdf5/serial/include
#/usr/include/hdf5/serial/
#/usr/lib/x86_64-linux-gnu/hdf5/serial
#/usr/lib/x86_64-linux-gnu/
#/usr/include/hdf5/serial/


RUN cd AxiSEM3D; cmake --build build;

ENV PATH="/home/axisem/AxiSEM3D/build:${PATH}"
