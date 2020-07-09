#
#  CMakeLists.txt
#  AxiSEM3D
#
#  Created by Kuangdai Leng on 8/15/19.
#  Copyright © 2019 Kuangdai Leng. All rights reserved.
#

#  CMakeLists.txt

########################################################################
###################### Build Settings (EDIT BELOW) #####################
########################################################################

# compiler suite
# ** either hardcode here or pass by -Dcc, -Dcxx and -Dftn after cmake
set(CMAKE_C_COMPILER            ${cc})
set(CMAKE_CXX_COMPILER          ${cxx})
set(CMAKE_Fortran_COMPILER      ${ftn})

# preferred installation prefix of external packages
# ** the roots can also be set by environment variables, e.g.,
#    $ export EIGEN3_ROOT=$HOME/Developer/eigen3
# ** by default, the following CMake variables are initialized by
#    the environment variables; if specified otherwise, they will
#    have priority over the environment variables
set(EIGEN3_ROOT                 $ENV{EIGEN3_ROOT})
set(BOOST_ROOT                  $ENV{BOOST_ROOT})
set(FFTW_ROOT                   $ENV{FFTW_ROOT})
set(METIS_ROOT                  $ENV{METIS_ROOT})
set(NETCDF_ROOT                 $ENV{NETCDF_ROOT})
# must use the same HDF5 used to build NetCDF
set(HDF5_ROOT                   $ENV{HDF5_ROOT})

# additional libraries to link to
# ** add libs if linking fails with missing symbols, e.g., "-lcurl"
set(ADDITIONAL_LIBS "-lcurl")

# polynomial order of spectral elements (from 1 to 8)
set(NPOL 4)

# solver precision
# ** using double precision increases memory usage but barely
#    affects solver performance
set(USE_DOUBLE false)

# memory saving mode
# ** minimize memory usage with some performance hit
set(MEMORY_SAVING_MODE false)

# serial build
# ** build AxiSEM3D without MPI
set(SERIAL_BUILD false)

# use parallel NetCDF for output
# ** parallel NetCDF is supported but not mandatory
set(USE_PARALLEL_NETCDF false)

########################################################################
###################### Build Settings (EDIT ABOVE) #####################
########################################################################



########### NOTE: Users are less likely to edit lines below ############



################# cmake setup #################
# use "Release" as the default build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# project
cmake_minimum_required(VERSION 3.15.3)
project(AxiSEM3D)
enable_language(Fortran)

# CMP0074 issue since CMake 3.12
if(POLICY CMP0074)
cmake_policy(SET CMP0074 NEW)
endif()


################# macros passed to AxiSEM3D #################
# version
add_definitions(-D_VERSION="1.0")

# source dir
add_definitions(-D_PROJ_DIR=\"${PROJECT_SOURCE_DIR}\")

# nPol
add_definitions(-D_NPOL=${NPOL})

# double
if(USE_DOUBLE)
add_definitions(-D_USE_DOUBLE)
endif()

# memory
if(MEMORY_SAVING_MODE)
add_definitions(-D_SAVE_MEMORY)
endif()

# serial build
if(SERIAL_BUILD)
add_definitions(-D_SERIAL_BUILD)
set(USE_PARALLEL_NETCDF false)
endif()

# parallel NetCDF
if(USE_PARALLEL_NETCDF)
add_definitions(-D_USE_PARALLEL_NETCDF)
endif()


################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# mpi
if(NOT SERIAL_BUILD)
find_package(MPI COMPONENTS C REQUIRED)
include_directories(${MPI_C_INCLUDE_DIRS})
endif()

# eigen
# FindEigen3.cmake only uses ENV in HINTS
set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
find_package(Eigen3 3.3.90 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# boost
# BOOST_ROOT can be case-sensitive depending on cmake version
set(Boost_ROOT ${BOOST_ROOT})
set(ENV{Boost_ROOT} $ENV{BOOST_ROOT})
find_package(Boost 1.71.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# fftw
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

# metis
# FindMETIS.cmake uses METIS_DIR instead of METIS_ROOT in HINTS
set(METIS_DIR ${METIS_ROOT})
set(ENV{METIS_DIR} $ENV{METIS_ROOT})
find_package(Metis REQUIRED QUIET)
include_directories(${METIS_INCLUDE_DIRS})

# NetCDF
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})
# check parallel
if(${USE_PARALLEL_NETCDF} AND (NOT EXISTS "${NETCDF_INCLUDES}/netcdf_par.h"))
message(FATAL_ERROR "Missing header netcdf_par.h under ${NETCDF_INCLUDES}/.\
\nUse a NetCDF build with parallel support or disable USE_PARALLEL_NETCDF.")
endif()

# hdf5
find_package(HDF5 COMPONENTS C COMPONENTS HL REQUIRED)


################# include paths #################
# find all directories containing .hpp
# 1) scan *.hpp files in src/
file(GLOB_RECURSE HPP_PATH_FILES src/*.hpp)
# 2) separate path and filename
set(HPP_PATHS "")
foreach(HPP_PATH_FILE ${HPP_PATH_FILES})
get_filename_component(HPP_PATH ${HPP_PATH_FILE} PATH)
set(HPP_PATHS ${HPP_PATHS} ${HPP_PATH})
endforeach()
# 3) remove duplicates
list(REMOVE_DUPLICATES HPP_PATHS)

# include all directories containing .hpp
include_directories(
# src
${HPP_PATHS}
# local externals
external
)


################# source files #################
# find *.cpp
file(GLOB_RECURSE CPP_FILES src/*.cpp)
# find *.c
file(GLOB_RECURSE C_FILES src/*.c)
# find *.f90
file(GLOB_RECURSE F90_FILES src/*.f90)

# add all source files
add_executable(
axisem3d
# src
${CPP_FILES}
${C_FILES}
${F90_FILES}
# non-header-only local externals
external/yaml/Yaml.cpp
)


################# linking #################
target_link_libraries(
axisem3d
${MPI_C_LIBRARIES}
${FFTW_LIBRARIES}
${METIS_LIBRARIES}
${NETCDF_LIBRARIES}
${HDF5_HL_LIBRARIES}
${HDF5_LIBRARIES}
${ADDITIONAL_LIBS}
)

# standard C++17
set_property(TARGET axisem3d PROPERTY CXX_STANDARD 17)


################# summary #################
# upper case of CMAKE_BUILD_TYPE for compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE_U)

# macro to message dependencies nicely
macro(msglist myhead mylist)
# length of myhead
string(LENGTH ${myhead} lenhead)
# indent
string(REPEAT " " ${lenhead} indent)
# replace ; with \n + indent
string(REPLACE ";" "\n${indent}" newlist "${mylist}")
message("${myhead}${newlist}")
endmacro()

# print summary
message(STATUS "")
message(STATUS "")
message(STATUS "")
message(STATUS "========================== Summary ===========================")
message("   AxiSEM3D options______________________________________________")
message("     SEM n-polynomial  =  ${NPOL}")
if(${USE_DOUBLE})
message("     Solver precision  =  double")
else()
message("     Solver precision  =  single")
endif()
message("     Memory-saving     =  ${MEMORY_SAVING_MODE}")
message("     Serial build      =  ${SERIAL_BUILD}")
message("     Parallel netcdf   =  ${USE_PARALLEL_NETCDF}")
message("   Build settings________________________________________________")
message("     Compilers")
message("       C               =  ${CMAKE_C_COMPILER}")
message("       C++             =  ${CMAKE_CXX_COMPILER}")
message("       Fortran         =  ${CMAKE_Fortran_COMPILER}")
message("     Compiler flags")
message("       C               =  ${CMAKE_C_FLAGS_${BUILD_TYPE_U}}")
message("       C++             =  ${CMAKE_CXX_FLAGS_${BUILD_TYPE_U}}")
message("       Fortran         =  ${CMAKE_Fortran_FLAGS_${BUILD_TYPE_U}}")
message("     CMake build type  =  ${CMAKE_BUILD_TYPE}")
message("     Source directory  =  ${PROJECT_SOURCE_DIR}")
message("   Dependencies__________________________________________________")
message("     * Eigen3")
msglist("       user-specified root: " "${EIGEN3_ROOT}")
msglist("       found include paths: " "${EIGEN3_INCLUDE_DIR}")
message("     * Boost")
msglist("       user-specified root: " "${BOOST_ROOT}")
msglist("       found include paths: " "${Boost_INCLUDE_DIR}")
message("     * FFTW3")
msglist("       user-specified root: " "${FFTW_ROOT}")
msglist("       found include paths: " "${FFTW_INCLUDES}")
msglist("       linked to libraries: " "${FFTW_LIBRARIES}")
message("     * METIS")
msglist("       user-specified root: " "${METIS_ROOT}")
msglist("       found include paths: " "${METIS_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${METIS_LIBRARIES}")
message("     * NetCDF")
msglist("       user-specified root: " "${NETCDF_ROOT}")
msglist("       found include paths: " "${NETCDF_INCLUDES}")
msglist("       linked to libraries: " "${NETCDF_LIBRARIES}")
message("     * HDF5")
msglist("       user-specified root: " "${HDF5_ROOT}")
msglist("       linked to libraries: " "${HDF5_LIBRARIES}")
msglist("                            " "${HDF5_HL_LIBRARIES}")
message("     * MPI")
if (${SERIAL_BUILD})
message("       This is a serial build without MPI.")
else()
if("${MPI_C_INCLUDE_DIRS}" STREQUAL "")
message("       Using MPI wrapper compilers, no external link required.")
else()
msglist("       found include paths: " "${MPI_C_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${MPI_C_LIBRARIES}")
endif()
endif()
message("     * Others")
message("       user-specified additional libraries: ${ADDITIONAL_LIBS}")
message(STATUS "==============================================================")
message(STATUS "")
message(STATUS "")
message(STATUS "")
