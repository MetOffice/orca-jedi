#
# (C) British Crown Copyright 2024 Met Office
#

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

list( APPEND oops_orcamodel_src_files
errorcovariance/ErrorCovariance.cc
errorcovariance/ErrorCovariance.h
geometry/GeometryParameters.h
geometry/Geometry.h
geometry/Geometry.cc
interpolator/InterpolatorParameters.h
interpolator/Interpolator.h
interpolator/Interpolator.cc
increment/Increment.cc
increment/Increment.h
model/ModelBiasCovariance.h
model/ModelBias.h
model/ModelBiasIncrement.h
model/Model.h
state/State.h
state/State.cc
state/StateIOUtils.h
state/StateIOUtils.cc
nemo_io/NemoFieldReader.h
nemo_io/NemoFieldReader.cc
nemo_io/NemoFieldWriter.h
nemo_io/NemoFieldWriter.cc
nemo_io/AtlasIndex.h
nemo_io/ReadServer.h
nemo_io/ReadServer.cc
nemo_io/WriteServer.h
nemo_io/WriteServer.cc
utilities/OrcaModelTraits.h
utilities/Types.h
utilities/Types.cc
variablechanges/VariableChangeParameters.h
variablechanges/VariableChange.h
)



ecbuild_add_library( TARGET   orcamodel
                     SOURCES  ${oops_orcamodel_src_files}
                     PUBLIC_LIBS  oops ufo atlas atlas-orca eckit netcdf
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${OOPS_LINKER_LANGUAGE}
                    )
                  #target_link_libraries( orcamodel PUBLIC MPI::MPI_Fortran MPI::MPI_CXX )
target_link_libraries( orcamodel PUBLIC NetCDF::NetCDF_CXX)
target_include_directories(orcamodel PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
