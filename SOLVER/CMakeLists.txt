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

# NOTE: these parameters can be either hardcoded in this box or
#       passed by -D (e.g., -Dcc=mpicc for CMAKE_C_COMPILER)


################## compiler options ##################
# compilers or compiler wrappers
set(CMAKE_C_COMPILER            ${cc})
set(CMAKE_CXX_COMPILER          ${cxx})
set(CMAKE_Fortran_COMPILER      ${ftn})

# compiler flags in addition to the standard ones (such as -O3 -DNDEBUG)
set(ADDITIONAL_COMPILER_FLAGS   ${flags})

# additional libraries to link to if linking fails with missing symbols
set(ADDITIONAL_LIBRARIES        ${libs})


################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN3_ROOT                 ${eigen})
set(BOOST_ROOT                  ${boost})
set(FFTW_ROOT                   ${fftw})
set(METIS_ROOT                  ${metis})
set(NETCDF_ROOT                 ${netcdf})
# omitting -Dhdf5 disables explicit linking to HDF5
set(HDF5_ROOT                   ${hdf5})

# use parallel NetCDF
# parallel NetCDF can enhance I/O performance and simplify
# post-processing, supported but not mandatory
set(USE_PARALLEL_NETCDF         ${par_netcdf})


################## solver options ##################
# polynomial order of spectral elements (from 1 to 8)
set(NPOL_SEM                    ${npol})

# solver precision
# using double precision increases memory usage but barely
# affects solver performance
set(USE_DOUBLE                  ${double})


########################################################################
###################### Build Settings (EDIT ABOVE) #####################
########################################################################



########### NOTE: Users are less likely to edit lines below ############



################# set default values #################
# macro to set default
macro(setDefault var val)
if(NOT ${var})
set(${var} ${val})
endif()
endmacro()

# set default values
setDefault(CMAKE_C_COMPILER            gcc)
setDefault(CMAKE_CXX_COMPILER          g++)
setDefault(CMAKE_Fortran_COMPILER      gfortran)
setDefault(ADDITIONAL_COMPILER_FLAGS   "")
setDefault(ADDITIONAL_LIBRARIES        "")
setDefault(EIGEN3_ROOT                 "<empty>")
setDefault(BOOST_ROOT                  "<empty>")
setDefault(FFTW_ROOT                   "<empty>")
setDefault(METIS_ROOT                  "<empty>")
setDefault(NETCDF_ROOT                 "<empty>")
setDefault(HDF5_ROOT                   "<empty>")
setDefault(USE_PARALLEL_NETCDF         false)
setDefault(NPOL_SEM                    4)
setDefault(USE_DOUBLE                  false)


################# cmake setup #################
# use "Release" as the default build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

# project
cmake_minimum_required(VERSION 3.10.0)
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
add_definitions(-D_NPOL=${NPOL_SEM})

# double
if(${USE_DOUBLE})
add_definitions(-D_USE_DOUBLE)
endif()

# memory
set(MEMORY_SAVING_MODE ${save_mem})
setDefault(MEMORY_SAVING_MODE false)
if(${MEMORY_SAVING_MODE})
add_definitions(-D_SAVE_MEMORY)
endif()

# serial build
set(SERIAL_BUILD ${serial})
setDefault(SERIAL_BUILD false)
if(${SERIAL_BUILD})
add_definitions(-D_SERIAL_BUILD)
set(USE_PARALLEL_NETCDF false)
endif()

# parallel NetCDF
if(${USE_PARALLEL_NETCDF})
add_definitions(-D_USE_PARALLEL_NETCDF)
endif()


################# find packages #################
# path of *.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# mpi
if(NOT ${SERIAL_BUILD})
find_package(MPI COMPONENTS C REQUIRED)
include_directories(${MPI_C_INCLUDE_DIRS})
endif()

# eigen
# FindEigen3.cmake only uses ENV in HINTS
set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
find_package(Eigen3 3.3.90 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# boost
set(ENV{BOOST_ROOT} ${BOOST_ROOT})
find_package(Boost 1.73.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# fftw
set(ENV{FFTW_ROOT} ${FFTW_ROOT})
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

# metis
# FindMETIS.cmake uses METIS_DIR instead of METIS_ROOT in HINTS
set(METIS_DIR ${METIS_ROOT})
set(ENV{METIS_DIR} ${METIS_ROOT})
set(ENV{METIS_ROOT} ${METIS_ROOT})
find_package(METIS REQUIRED)
include_directories(${METIS_INCLUDE_DIRS})

# NetCDF
set(ENV{NETCDF_ROOT} ${NETCDF_ROOT})
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})
# check parallel
if(${USE_PARALLEL_NETCDF})
if(NOT EXISTS "${NETCDF_INCLUDES}/netcdf_par.h")
message(FATAL_ERROR "\nMissing header netcdf_par.h under ${NETCDF_INCLUDES}/.\
\nUse a NetCDF build with parallel support or disable parallel NetCDF by \
-Dpar_netcdf=false.")
endif()
endif()

# hdf5
if(hdf5)
set(ENV{HDF5_ROOT} ${HDF5_ROOT})
find_package(HDF5 COMPONENTS C COMPONENTS HL REQUIRED)
endif()


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
# compiler flags
# upper case of CMAKE_BUILD_TYPE for compiler flags
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD)
set(ADD_FLAGS ${ADDITIONAL_COMPILER_FLAGS})
set(CMAKE_C_FLAGS_${BUILD} "${CMAKE_C_FLAGS_${BUILD}} ${ADD_FLAGS}")
set(CMAKE_CXX_FLAGS_${BUILD} "${CMAKE_CXX_FLAGS_${BUILD}} ${ADD_FLAGS}")
set(CMAKE_Fortran_FLAGS_${BUILD} "${CMAKE_Fortran_FLAGS_${BUILD}} ${ADD_FLAGS}")

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
${ADDITIONAL_LIBRARIES}
)

# standard C++17
set_property(TARGET axisem3d PROPERTY CXX_STANDARD 17)


################# summary #################
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
message("   Build settings________________________________________________")
message("     Compilers")
message("       C               =  ${CMAKE_C_COMPILER}")
message("       C++             =  ${CMAKE_CXX_COMPILER}")
message("       Fortran         =  ${CMAKE_Fortran_COMPILER}")
message("     Compiler flags")
message("       C               =  ${CMAKE_C_FLAGS_${BUILD}}")
message("       C++             =  ${CMAKE_CXX_FLAGS_${BUILD}}")
message("       Fortran         =  ${CMAKE_Fortran_FLAGS_${BUILD}}")
message("     CMake build type  =  ${CMAKE_BUILD_TYPE}")
message("     Source directory  =  ${PROJECT_SOURCE_DIR}")
message("   Dependencies__________________________________________________")
message("     * MPI")
if(${SERIAL_BUILD})
message("       This is a serial build without MPI.")
else()
if("${MPI_C_INCLUDE_DIRS}" STREQUAL "")
message("       Using MPI wrapper compilers, no explicit link required.")
else()
msglist("       found include paths: " "${MPI_C_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${MPI_C_LIBRARIES}")
endif()
endif()
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
msglist("       use parallel NetCDF: " "${USE_PARALLEL_NETCDF}")
message("     * HDF5")
if(NOT hdf5)
message("       Linking to HDF5 is disabled by omitting -Dhdf5.")
else()
if("${HDF5_LIBRARIES}" STREQUAL "")
message("       Using HDF5 wrapper compilers, no explicit link required.")
else()
msglist("       user-specified root: " "${HDF5_ROOT}")
msglist("       linked to libraries: " "${HDF5_LIBRARIES}")
msglist("                            " "${HDF5_HL_LIBRARIES}")
endif()
endif()
message("     * Others")
if("${ADDITIONAL_LIBRARIES}" STREQUAL "")
message("       No user-specified additional libraries.")
else()
message("       user-specified additional libraries: ${ADDITIONAL_LIBRARIES}")
endif()
message("   AxiSEM3D options______________________________________________")
message("     SEM n-polynomial  =  ${NPOL_SEM}")
if(${USE_DOUBLE})
message("     Solver precision  =  double")
else()
message("     Solver precision  =  single")
endif()
message(STATUS "==============================================================")
message(STATUS "")
message(STATUS "")
message(STATUS "")
