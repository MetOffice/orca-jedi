include_directories( ${PROJECT_SOURCE_DIR}/src )

#target_link_libraries( meara_test_nemo_feedback NetCDF::NetCDF_CXX )

ecbuild_add_test( TARGET  test_orcajedi_nemo_io_field_reader.x
                  SOURCES test_nemo_io_field_reader.cc
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_nemo_io_field_writer.x
                  SOURCES test_nemo_io_field_writer.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_nemo_io_orca_field_reader.x
                  SOURCES test_nemo_io_orca_field_reader.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_nemo_io_read_server.x
                  SOURCES test_nemo_io_read_server.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_geometry.x
                  SOURCES test_geometry.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_state.x
                  SOURCES test_state.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_interpolator.x
                  SOURCES test_interpolator.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  test_orcajedi_increment.x
                  SOURCES test_increment.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )
