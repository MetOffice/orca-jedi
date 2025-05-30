# (C) British Crown Copyright 2024 Met Office

############################################################################################
# orca-jedi

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )

################################################################################
# Initialise project orca-jedi

project( orca-jedi LANGUAGES CXX Fortran )
ecbuild_declare_project()

set( CMAKE_DIRECTORY_LABELS "orca-jedi" )

################################################################################
# Required packages

## Dependencies
find_package( jedicmake QUIET )  # Prefer find modules from jedi-cmake

find_package( NetCDF COMPONENTS CXX)
ecbuild_debug( "   NetCDF_FEATURES: [${NetCDF_FEATURES}]" )

find_package( eckit 1.18 COMPONENTS LZ4 REQUIRED )
ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

find_package( atlas 0.38.1 REQUIRED )
ecbuild_debug( "   atlas_FEATURES: [${atlas_FEATURES}]" )

find_package( atlas-orca REQUIRED )
ecbuild_debug( "   atlas_orca_FEATURES: [${atlas_orca_FEATURES}]" )

find_package( oops REQUIRED )
ecbuild_debug( "   oops_FEATURES : [${oops_FEATURES}]" )

find_package( saber REQUIRED)
ecbuild_debug( "   saber_FEATURES : [${saber_FEATURES}]" )

find_package( ufo REQUIRED)
ecbuild_debug( "   ufo_FEATURES : [${ufo_FEATURES}]" )

find_package( nemo-feedback )

################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

# check if Doxygen is installed
find_package(Doxygen)
ecbuild_add_option( FEATURE ORCA_JEDI_DOCS
                    DESCRIPTION "orca-jedi documentation"
                    DEFAULT OFF
                    REQUIRED_PACKAGES "Doxygen" )

################################################################################
# sources

#include( atlas_compile_flags )

include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${CMAKE_CURRENT_BINARY_DIR}/src
                    )

add_subdirectory( src )

add_subdirectory ( docs )

################################################################################
# Export and summarize

ecbuild_install_project( NAME orcamodel )

ecbuild_print_summary()

