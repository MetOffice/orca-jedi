include_directories( ${PROJECT_SOURCE_DIR}/src )

ecbuild_add_test( TARGET  orcajedi_nemo_io_field_reader.x
                  SOURCES test_nemo_io_field_reader.cc
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_nemo_io_field_writer.x
                  SOURCES test_nemo_io_field_writer.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_nemo_io_read_server_MPI2.x
                  SOURCES test_nemo_io_read_server.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_nemo_io_write_server_MPI2.x
                  SOURCES test_nemo_io_write_server.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  MPI        2
                  CONDITION  eckit_HAVE_MPI
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_nemo_io_write_server.x
                  SOURCES test_nemo_io_write_server.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_geometry.x
                  SOURCES test_geometry.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_state.x
                  SOURCES test_state.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  OMP     4
                  CONDITION  eckit_HAVE_OMP
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_interpolator.x
                  SOURCES test_interpolator.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_increment.x
                  SOURCES test_increment.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_regridder.x
                  SOURCES test_regridder.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )

ecbuild_add_test( TARGET  orcajedi_source_extender.x
                  SOURCES test_source_extender.cc
                  ENVIRONMENT ${ATLAS_TEST_ENVIRONMENT}
                  LIBS    orcamodel )
