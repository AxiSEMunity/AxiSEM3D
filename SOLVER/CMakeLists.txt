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
#       passed by -D (e.g., -Dcxx=mpicxx for CMAKE_CXX_COMPILER)


################## compiler options ##################
# compiler
if(DEFINED cxx AND NOT cxx STREQUAL "")
  set(CMAKE_CXX_COMPILER "${cxx}" CACHE FILEPATH "C++ compiler" FORCE)
endif()

# additional compiler flags
set(ADDITIONAL_CXX_FLAGS "${flags}")

# additional link options
set(ADDITIONAL_LINK_OPTIONS "${links}")


################## dependencies ##################
# preferred installation prefix of dependencies
set(EIGEN3_ROOT                 ${eigen})
set(BOOST_ROOT                  ${boost})
set(FFTW_ROOT                   ${fftw})
set(METIS_ROOT                  ${metis})
set(NETCDF_ROOT                 ${netcdf})
set(HDF5_ROOT                   ${hdf5})

# use parallel NetCDF
# parallel NetCDF is supported but not mandatory
set(USE_PARALLEL_NETCDF         ${par_netcdf})


################## solver options ##################
# polynomial order of spectral elements (from 1 to 8)
set(NPOL_SEM                    ${npol})

# solver precision
# using double precision doubles runtime memory and output size but
# affects performance only slightly on a modern 64-bit architecture
set(USE_DOUBLE                  ${double})


########################################################################
###################### Build Settings (EDIT ABOVE) #####################
########################################################################



########### NOTE: Users are less likely to edit lines below ############



################# set default values #################
# macro to set default
macro(setDefault var val)
  if(NOT DEFINED "${var}")
    set("${var}" "${val}")
  endif()
endmacro()

# set default values
setDefault(ADDITIONAL_CXX_FLAGS        "")
setDefault(ADDITIONAL_LINK_OPTIONS     "")
setDefault(EIGEN3_ROOT                 "<empty>")
setDefault(BOOST_ROOT                  "<empty>")
setDefault(FFTW_ROOT                   "<empty>")
setDefault(METIS_ROOT                  "<empty>")
setDefault(NETCDF_ROOT                 "<empty>")
setDefault(HDF5_ROOT                   "<empty>")
setDefault(USE_PARALLEL_NETCDF         OFF)
setDefault(NPOL_SEM                    4)
setDefault(USE_DOUBLE                  OFF)


################# cmake setup #################
# project
cmake_minimum_required(VERSION 4.2)
project(AxiSEM3D)

# use "Release" as the default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()


################# macros passed to AxiSEM3D #################
# version
add_definitions(-D_VERSION="1.0")

# source dir
add_definitions(-D_PROJ_DIR=\"${PROJECT_SOURCE_DIR}\")

# nPol
add_definitions(-D_NPOL=${NPOL_SEM})

# double
if(USE_DOUBLE)
  add_definitions(-D_USE_DOUBLE)
endif()

# memory
set(MEMORY_SAVING_MODE "${save_mem}")
setDefault(MEMORY_SAVING_MODE OFF)
if(MEMORY_SAVING_MODE)
  add_definitions(-D_SAVE_MEMORY)
endif()

# serial build
set(SERIAL_BUILD "${serial}")
setDefault(SERIAL_BUILD OFF)
if(SERIAL_BUILD)
  add_definitions(-D_SERIAL_BUILD)
  set(USE_PARALLEL_NETCDF OFF)
endif()

# parallel NetCDF
if(USE_PARALLEL_NETCDF)
  add_definitions(-D_USE_PARALLEL_NETCDF)
endif()

# skip _MM_SET_FLUSH_ZERO_MODE in fenv.cpp
set(SKIP_MM_SET "${skip_mm_set}")
setDefault(SKIP_MM_SET OFF)
if(SKIP_MM_SET)
  add_definitions(-D_SKIP_DISABLE_SSE_DENORMS)
endif()


################# find packages #################
# path of *.cmake files
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# mpi
if(NOT SERIAL_BUILD)
  find_package(MPI COMPONENTS CXX REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_DIRS})
endif()

# eigen
find_package(Eigen3 3.4 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# boost
set(BOOST_INCLUDEDIR "${BOOST_ROOT}/include")
set(BOOST_LIBRARYDIR "${BOOST_ROOT}/lib")
find_package(Boost 1.85 REQUIRED)

# fftw
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDES})

# metis
# FindMETIS.cmake uses METIS_DIR instead of METIS_ROOT in HINTS
set(METIS_DIR "${METIS_ROOT}")
find_package(METIS REQUIRED)
include_directories(${METIS_INCLUDE_DIRS})

# NetCDF
find_package(NetCDF REQUIRED)
include_directories(${NETCDF_INCLUDES})
# check parallel
if(USE_PARALLEL_NETCDF)
  if(NOT EXISTS "${NETCDF_INCLUDES}/netcdf_par.h")
    message(FATAL_ERROR "\nMissing header netcdf_par.h under ${NETCDF_INCLUDES}/.\
\nUse a NetCDF build with parallel support or disable parallel NetCDF by \
-Dpar_netcdf=false.")
  endif()
endif()

# hdf5
if(NOT HDF5_ROOT STREQUAL "")
  find_package(HDF5 REQUIRED COMPONENTS C HL)
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
# find *.cpp
file(GLOB_RECURSE CPP_FILES src/*.cpp)

# add all source files
add_executable(
axisem3d
# src
${CPP_FILES}
# non-header-only local externals
external/yaml/Yaml.cpp
)

# standard C++17
set_property(TARGET axisem3d PROPERTY CXX_STANDARD 17)

# compiler flags
# user-provided extra flags (space-separated string supported)
if(NOT "${ADDITIONAL_CXX_FLAGS}" STREQUAL "")
  separate_arguments(_A3D_USER_CXX_FLAGS NATIVE_COMMAND "${ADDITIONAL_CXX_FLAGS}")
  target_compile_options(axisem3d PRIVATE ${_A3D_USER_CXX_FLAGS})
endif()

# suppress #pragma warnings
target_compile_options(axisem3d PRIVATE -Wno-unknown-pragmas)


################# linking #################
if(NOT HDF5_ROOT STREQUAL "")
  target_link_libraries(
    axisem3d
    ${MPI_CXX_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${METIS_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${HDF5_LIBRARIES}
  )
else()
  target_link_libraries(
    axisem3d
    ${MPI_CXX_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${METIS_LIBRARIES}
    ${NETCDF_LIBRARIES}
  )
endif()

# additional user link options (space-separated string supported)
if(NOT "${ADDITIONAL_LINK_OPTIONS}" STREQUAL "")
  separate_arguments(_A3D_USER_LINK_OPTIONS NATIVE_COMMAND "${ADDITIONAL_LINK_OPTIONS}")
  target_link_options(axisem3d PRIVATE ${_A3D_USER_LINK_OPTIONS})
endif()


################# summary #################
# macro to message dependencies nicely
macro(msglist myhead mylist)
# form indent string
string(LENGTH ${myhead} lenhead)
# this neat solution requires 3.15
# string(REPEAT " " ${lenhead} indent)
set(indent "")
foreach(idx RANGE 1 ${lenhead})
set(indent "${indent} ")
endforeach()
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
message("     C++ compiler      =  ${CMAKE_CXX_COMPILER}")
message("     Compiler flags    =  ${CMAKE_CXX_FLAGS_${BUILD}}")
message("     Product type      =  ${CMAKE_BUILD_TYPE}")
message("     Source directory  =  ${PROJECT_SOURCE_DIR}")
message("   Dependencies__________________________________________________")
message("     * MPI")
if(${SERIAL_BUILD})
message("       This is a serial build without MPI.")
else()
if("${MPI_CXX_INCLUDE_DIRS}" STREQUAL "")
message("       Using MPI wrapper compilers, no explicit link required.")
else()
msglist("       found include paths: " "${MPI_CXX_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${MPI_CXX_LIBRARIES}")
endif()
endif()
message("     * Eigen3")
msglist("       user-specified root: " "${EIGEN3_ROOT}")
msglist("       found include paths: " "${EIGEN3_INCLUDE_DIR}")
message("     * Boost")
msglist("       user-specified root: " "${BOOST_ROOT}")
msglist("       found include paths: " "${Boost_INCLUDE_DIRS}")
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
if(NOT DEFINED hdf5)
message("       Linking to HDF5 is disabled by omitting -Dhdf5.")
else()
msglist("       user-specified root: " "${HDF5_ROOT}")
if("${HDF5_LIBRARIES}" STREQUAL "")
message("       Using HDF5 wrapper compilers, no explicit link required.")
else()
msglist("       linked to libraries: " "${HDF5_HL_LIBRARIES}")
msglist("                            " "${HDF5_LIBRARIES}")
endif()
endif()
message("     * Link options")
if("${LINK_OPTIONS}" STREQUAL "")
message("       No user-specified link options.")
else()
message("       user-specified link options: ${LINK_OPTIONS}")
endif()
message("   Solver options________________________________________________")
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
