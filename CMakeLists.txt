# SPDX-License-Identifier: MIT
# Copyright (C) 2019 - 2026 by the AxiSEM3D authors
#
# This file is part of the AxiSEM3D library. See the LICENSE file for details.

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
# parallel NetCDF is supported but not mandatory
set(USE_PARALLEL_NETCDF OFF CACHE BOOL "If ON, use parallel NetCDF. Optional.")

################## solver options ##################
set(NPOL_SEM "4" CACHE STRING "polynomial order of spectral elements (from 1 to 8)")

# solver precision
# using double precision doubles runtime memory and output size but
# affects performance only slightly on a modern 64-bit architecture
set(USE_DOUBLE OFF CACHE BOOL "If ON use double precision, otherwise single")


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


################# cmake setup #################
# project
cmake_minimum_required(VERSION 3.26)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not allowed. Use a separate build directory, for example:\n"
    "  cmake -S ${CMAKE_SOURCE_DIR} -B build\n"
    "  cmake --build build")
endif()

project(AxiSEM3D)
include(CTest)

# use "Release" as the default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()


################# macros passed to AxiSEM3D #################
# grab version from VERSION file:
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" AXISEM3D_VERSION LIMIT_COUNT 1)
add_definitions(-DAXISEM3D_VERSION="${AXISEM3D_VERSION}")

# source dir
add_definitions(-DAXISEM3D_PROJ_DIR=\"${PROJECT_SOURCE_DIR}\")

# nPol
add_definitions(-D AXISEM3D_NPOL=${NPOL_SEM})

# double
if(USE_DOUBLE)
  add_definitions(-DAXISEM3D_USE_DOUBLE)
endif()

# memory
set(MEMORY_SAVING_MODE OFF CACHE BOOL "If ON, try to save system memory")
if(MEMORY_SAVING_MODE)
  add_definitions(-DAXISEM3D_SAVE_MEMORY)
endif()

# serial build
set(SERIAL_BUILD OFF CACHE BOOL "If ON, perform a serial build without MPI.")
if(SERIAL_BUILD)
  add_definitions(-DAXISEM3D_SERIAL_BUILD)
  set(USE_PARALLEL_NETCDF OFF CACHE FORCE)
endif()

# parallel NetCDF
if(USE_PARALLEL_NETCDF)
  add_definitions(-DAXISEM3D_USE_PARALLEL_NETCDF)
endif()

# skip _MM_SET_FLUSH_ZERO_MODE in fenv.cpp
set(SKIP_MM_SET OFF CACHE BOOL "if ON, skip _MM_SET_FLUSH_ZERO_MODE in fenv.cpp")
if(SKIP_MM_SET)
  add_definitions(-DAXISEM3D_SKIP_DISABLE_SSE_DENORMS)
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
set(EIGEN3_DIR "$ENV{EIGEN3_DIR}" CACHE PATH "An optional hint to a Eigen3 installation")
find_package(Eigen3 3.4 HINTS ${EIGEN3_DIR})
if (NOT Eigen3_FOUND)
  message(STATUS "Eigen3 CMake package not found; using FindEigen3.cmake fallback...")
  find_package(Eigen3 REQUIRED MODULE)
endif()
include_directories(${EIGEN3_INCLUDE_DIR})

# boost
# Silence a warning by CMake, the new policy is fine for
# us as boost is new enough to provide the necessary info.
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()
set(BOOST_DIR "$ENV{BOOST_DIR}" CACHE PATH "An optional hint to a Boost installation")
find_package(Boost 1.83 REQUIRED HINTS ${BOOST_DIR})
include_directories(${Boost_INCLUDE_DIRS})

# FFTW3 (optional prefix: -DFFTW3_DIR= or env FFTW3_DIR; legacy -Dfftw= still accepted)
set(FFTW3_DIR "$ENV{FFTW3_DIR}" CACHE PATH "Optional hint to an FFTW3 installation prefix")
if(DEFINED fftw AND NOT "${fftw}" STREQUAL "")
  set(FFTW3_DIR "${fftw}" CACHE PATH "Optional hint to an FFTW3 installation prefix" FORCE)
endif()
find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})

# metis
set(METIS_DIR "$ENV{METIS_DIR}" CACHE PATH "An optional hint to a METIS installation")

# metis: try CMake package config (some CMake-built installs); fall back to FindMETIS.cmake
# User prefix is METIS_ROOT (e.g. -Dmetis=$CONDA_PREFIX). Do not pre-set METIS_DIR: CMake sets
# METIS_DIR from find_package(METIS CONFIG) and may leave METIS_DIR-NOTFOUND on failure, which
# would break FindMETIS.cmake unless we override before the MODULE search.
set(METIS_DIR_BACKUP "${METIS_DIR}")
find_package(METIS QUIET CONFIG HINTS "${METIS_DIR}" PATHS "${METIS_DIR}" "$ENV{METIS_DIR}")
if(NOT METIS_FOUND)
  message(STATUS "METIS CMake package not found; using FindMETIS.cmake fallback...")
  set(METIS_DIR "${METIS_DIR_BACKUP}" CACHE PATH "" FORCE)
  # MODULE mode does not accept HINTS/PATHS (CMake 4+); hints are via METIS_DIR cache above.
  find_package(METIS REQUIRED MODULE)
endif()
# Map config-package names to METIS_* used below (FindMETIS.cmake sets these directly).
if(NOT METIS_INCLUDE_DIRS)
  if(METIS_INCLUDE_DIR)
    set(METIS_INCLUDE_DIRS "${METIS_INCLUDE_DIR}")
  elseif(TARGET METIS::metis)
    get_target_property(_metis_iface_includes METIS::metis INTERFACE_INCLUDE_DIRECTORIES)
    if(_metis_iface_includes AND NOT _metis_iface_includes MATCHES "-NOTFOUND$")
      set(METIS_INCLUDE_DIRS "${_metis_iface_includes}")
    endif()
  endif()
endif()
if(NOT METIS_LIBRARIES)
  if(TARGET METIS::metis)
    set(METIS_LIBRARIES METIS::metis)
  elseif(METIS_LIBRARY)
    set(METIS_LIBRARIES "${METIS_LIBRARY}")
  endif()
endif()
include_directories(${METIS_INCLUDE_DIRS})

# NetCDF: try CMake package config first (conda, CMake-built installs); fall back to FindNetCDF.cmake
# for typical autotools installs without netCDFConfig.cmake.
set(NETCDF_DIR "" CACHE PATH "An optional hint to a NetCDF installation")
find_package(netCDF QUIET CONFIG HINTS "$ENV{NETCDF_DIR}" "${NETCDF_DIR}" PATHS "${NETCDF_DIR}")
if(NOT netCDF_FOUND)
  message(STATUS "NetCDF CMake package not found; using FindNetCDF.cmake fallback...")
  find_package(NetCDF REQUIRED MODULE)
endif()
# Map config-package names to NETCDF_* used below (FindNetCDF.cmake sets NETCDF_* directly).
if(NOT NETCDF_INCLUDES AND netCDF_INCLUDE_DIR)
  set(NETCDF_INCLUDES "${netCDF_INCLUDE_DIR}")
endif()
if(NOT NETCDF_LIBRARIES AND netCDF_LIBRARIES)
  set(NETCDF_LIBRARIES "${netCDF_LIBRARIES}")
endif()
include_directories(${NETCDF_INCLUDES})

# check parallel
if(USE_PARALLEL_NETCDF)
  if(NOT EXISTS "${NETCDF_INCLUDES}/netcdf_par.h")
    message(FATAL_ERROR "\nMissing header netcdf_par.h under ${NETCDF_INCLUDES}/.\
\nUse a NetCDF build with parallel support or disable parallel NetCDF by \
-DUSE_PARALLEL_NETCDF=OFF.")
  endif()
endif()

# hdf5
set(WITH_HDF5 ON CACHE BOOL "Build with HDF5 support")
if(WITH_HDF5)
  find_package(HDF5 COMPONENTS C HL)
  if(NOT HDF5_FOUND)
    message(WARNING "HDF5 not found, building without HDF5 support.")
    set(WITH_HDF5 OFF CACHE FORCE)
  else()
    include_directories(${HDF5_INCLUDE_DIRS})
    link_directories(${HDF5_LIBRARY_DIRS})
  endif()
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

# add all source files to the debug executable
add_executable(
axisem3d-debug
# src
${CPP_FILES}
# non-header-only local externals
external/yaml/Yaml.cpp
)

# standard C++17
set_property(TARGET axisem3d PROPERTY CXX_STANDARD 17)
set_property(TARGET axisem3d-debug PROPERTY CXX_STANDARD 17)


# compiler flags
# user-provided extra flags (space-separated string supported)
if(NOT "${ADDITIONAL_CXX_FLAGS}" STREQUAL "")
  separate_arguments(_A3D_USER_CXX_FLAGS NATIVE_COMMAND "${ADDITIONAL_CXX_FLAGS}")
  target_compile_options(axisem3d PRIVATE ${_A3D_USER_CXX_FLAGS})
  target_compile_options(axisem3d-debug PRIVATE ${_A3D_USER_CXX_FLAGS})
endif()

# Make sure to compile the debug target with debug symbols
# and no optimizations
target_compile_options(axisem3d-debug PRIVATE -g -O0)

# suppress #pragma warnings
target_compile_options(axisem3d PRIVATE -Wno-unknown-pragmas)
target_compile_options(axisem3d-debug PRIVATE -Wno-unknown-pragmas)


################# linking #################
target_link_libraries(
    axisem3d
    ${MPI_CXX_LIBRARIES}
    ${FFTW3_LIBRARIES}
    ${METIS_LIBRARIES}
    ${NETCDF_LIBRARIES}
  )

target_link_libraries(
    axisem3d-debug
    ${MPI_CXX_LIBRARIES}
    ${FFTW3_LIBRARIES}
    ${METIS_LIBRARIES}
    ${NETCDF_LIBRARIES}
  )

if(WITH_HDF5)
  target_link_libraries(axisem3d ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
  target_link_libraries(axisem3d-debug ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
endif()

# additional user link options (space-separated string supported)
if(NOT "${ADDITIONAL_LINK_OPTIONS}" STREQUAL "")
  separate_arguments(_A3D_USER_LINK_OPTIONS NATIVE_COMMAND "${ADDITIONAL_LINK_OPTIONS}")
  target_link_options(axisem3d PRIVATE ${_A3D_USER_LINK_OPTIONS})
  target_link_options(axisem3d-debug PRIVATE ${_A3D_USER_LINK_OPTIONS})
endif()

# target to run tools/indent.sh (clang-format) from project root
add_custom_target(indent
  COMMAND sh "${CMAKE_SOURCE_DIR}/tools/indent.sh"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMENT "Indenting source with clang-format (tools/indent.sh)"
)


# Integration tests
find_program(NUMDIFF_EXECUTABLE numdiff QUIET)
if(NOT NUMDIFF_EXECUTABLE)
  message(STATUS "numdiff not found, skipping integration tests. See README.md for installation instructions.")
endif()

if(NUMDIFF_EXECUTABLE AND EXISTS "${PROJECT_SOURCE_DIR}/tests")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)

  message(STATUS "Searching for tests in ${PROJECT_SOURCE_DIR}/tests:")
  file(GLOB _tests RELATIVE "" "${PROJECT_SOURCE_DIR}/tests/*/config.txt")

  set(_n_tests "0")
  list(SORT _tests)

  foreach(_test ${_tests})
    string(REPLACE "/config.txt" "" _testfolder ${_test})
    get_filename_component(_testname ${_testfolder} NAME)

    message(STATUS "  Test ${_testname}:")
    message(STATUS "    Path ${_testfolder}")

    if(NOT IS_DIRECTORY "${_testfolder}/input" OR
      NOT EXISTS "${_testfolder}/screen-output")
      message(STATUS "    Skipping integration test case '${_testname}': missing input/ or screen-output")
      continue()
    endif()

    # grab MPI ranks from config.txt if present:
    file(STRINGS "${_testfolder}/config.txt" _n_mpi_ranks REGEX "n_mpi_ranks[ ]+[0-9]+")
    if(_n_mpi_ranks)
      list(GET _n_mpi_ranks 0 _line)
      string(REGEX REPLACE "n_mpi_ranks[ ]+([0-9]+).*" "\\1" _n_mpi_ranks "${_line}")
    else()
      set(_n_mpi_ranks "1")
    endif()

    message(STATUS "    Running with ${_n_mpi_ranks} MPI ranks")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname})

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output
      COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/runs
      COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/output
      COMMAND echo mpirun -np ${_n_mpi_ranks} ${CMAKE_BINARY_DIR}/axisem3d --input ${_testfolder}/input --output ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/output
      COMMAND OMP_NUM_THREADS=1 mpirun -np ${_n_mpi_ranks} ${CMAKE_BINARY_DIR}/axisem3d --input ${_testfolder}/input --output ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/output 2>&1 | tee ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/raw-screen-output
      COMMAND ${CMAKE_SOURCE_DIR}/cmake/normalize_screen_output.sh ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/raw-screen-output ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output
      COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/runs
      DEPENDS axisem3d ${_test} ${CMAKE_SOURCE_DIR}/CMakeLists.txt
    )

    add_custom_target(tests.${_testname}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output)

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output.diff
    COMMAND test -f ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/runs
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/diff_test.sh
            ${NUMDIFF_EXECUTABLE}
            ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output
            ${_testfolder}/screen-output
    DEPENDS ${_testfolder}/screen-output ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output)
    add_custom_target(tests.${_testname}.screen-output.diff
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output.diff)

    add_dependencies(tests.${_testname} tests.${_testname}.screen-output.diff)

    # find all files under reference/; paths are relative to reference/
    file(GLOB_RECURSE _reference_files
         RELATIVE "${_testfolder}/reference"
         "${_testfolder}/reference/*")
    # For each reference file, create a comand to compare against output file
    foreach(_reference_file ${_reference_files})
      message(STATUS "    Reference file: ${_reference_file}")
      # name without slashes:
      string(REPLACE "/" "_" _nice_name ${_reference_file})

      add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/${_nice_name}.diff
      COMMAND ${CMAKE_SOURCE_DIR}/cmake/diff_test.sh
              ${NUMDIFF_EXECUTABLE}
              ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/output/${_reference_file}
              ${_testfolder}/reference/${_reference_file}
      DEPENDS ${_testfolder}/reference/${_reference_file} ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/screen-output)
      add_custom_target(tests.${_testname}.${_nice_name}.diff
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/tests/${_testname}/${_nice_name}.diff)
      add_dependencies(tests.${_testname} tests.${_testname}.${_nice_name}.diff)
    endforeach()

    add_test(NAME ${_testname}
    COMMAND ${CMAKE_COMMAND}
            -DBINARY_DIR=${CMAKE_BINARY_DIR}
            -DTESTNAME=tests.${_testname}
            -DERROR="Test ${_testname} failed"
            -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(${_testname} PROPERTIES TIMEOUT 600)

    math(EXPR _n_tests "${_n_tests}+1")
  endforeach()

  message(STATUS "Found ${_n_tests} tests.")
  message(STATUS "")
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
message("   AxiSEM3D version________________________________________________")
message("     Version           =  ${AXISEM3D_VERSION}")
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
msglist("       user-specified root: " "${EIGEN3_DIR}")
msglist("       found include paths: " "${EIGEN3_INCLUDE_DIR}")
message("     * Boost")
msglist("       user-specified root: " "${BOOST_DIR}")
msglist("       found include paths: " "${Boost_INCLUDE_DIRS}")
message("     * FFTW3")
msglist("       user-specified root: " "${FFTW3_DIR}")
msglist("       found include paths: " "${FFTW3_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${FFTW3_LIBRARIES}")
message("     * METIS")
msglist("       user-specified root: " "${METIS_DIR}")
msglist("       found include paths: " "${METIS_INCLUDE_DIRS}")
msglist("       linked to libraries: " "${METIS_LIBRARIES}")
message("     * NetCDF")
msglist("       user-specified root: " "${NETCDF_DIR}")
msglist("       found include paths: " "${NETCDF_INCLUDES}")
msglist("       linked to libraries: " "${NETCDF_LIBRARIES}")
msglist("       use parallel NetCDF: " "${USE_PARALLEL_NETCDF}")
message("     * HDF5 (${WITH_HDF5})")
if(NOT WITH_HDF5)
message("       Linking to HDF5 is disabled.")
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
