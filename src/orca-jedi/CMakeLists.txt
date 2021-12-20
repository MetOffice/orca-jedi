#
# (C) British Crown Copyright 2020 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_orcamodel_src_files
errorcovariance/ErrorCovariance.cc
errorcovariance/ErrorCovariance.h
geometry/GeometryParameters.h
geometry/Geometry.h
geometry/Geometry.cc
getvalues/GetValuesParameters.h
getvalues/GetValues.h
getvalues/GetValues.cc
getvalues/LinearGetValues.h
increment/Increment.cc
increment/Increment.h
model/ModelBiasCovariance.h
model/ModelBias.h
model/ModelBiasIncrement.h
model/Model.h
state/State.h
state/State.cc
state/StateReadUtils.h
state/StateReadUtils.cc
nemo_io/NemoFieldReader.h
nemo_io/NemoFieldReader.cc
utilities/OrcaModelTraits.h
variablechanges/VariableChangeParameters.h
variablechanges/LinearVariableChangeParameters.h
variablechanges/VariableChange.h
variablechanges/LinearVariableChange.h
)



ecbuild_add_library( TARGET   orcamodel
                     SOURCES  ${oops_orcamodel_src_files}
                     PUBLIC_LIBS  oops saber ufo atlas netcdf 
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
                  #target_link_libraries( orcamodel PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( orcamodel PUBLIC NetCDF::NetCDF_CXX)
target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
