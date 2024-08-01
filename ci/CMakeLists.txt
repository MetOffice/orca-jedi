# (C) British Crown Copyright 2024 Met Office

cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

project( orca-jedi-ci VERSION 1.0 LANGUAGES C CXX Fortran )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )

LIST( APPEND ATLAS_TEST_ENVIRONMENT
    ATLAS_ORCA_CACHING=1
    ATLAS_CACHE_PATH=${CMAKE_BINARY_DIR}/share
    PLUGINS_MANIFEST_PATH=${CMAKE_BINARY_DIR}/share/plugins )

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jedicmake" EXCLUDE_FROM_ALL)
if(NOT DEFINED jedicmake_DIR)
  set(jedicmake_DIR "${CMAKE_BINARY_DIR}/jedicmake")
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/oops" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ioda" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ufo" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/orca-jedi")
enable_testing()
