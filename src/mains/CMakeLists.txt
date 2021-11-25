if( ${nemo-feedback_FOUND} )
  message("      nemo-feedback_FOUND adding NEMO Feedback Writer")
  ecbuild_add_executable( TARGET  orcamodel_hofx.x
                          SOURCES orcamodelHofX.cc
                          LIBS    orcamodel nemo_feedback
                          DEFINITIONS NEMO_FEEDBACK_EXISTS
                        )

  ecbuild_add_executable( TARGET  orcamodel_hofx3D.x
                          SOURCES orcamodelHofX3D.cc
                          LIBS    orcamodel nemo_feedback
                          DEFINITIONS NEMO_FEEDBACK_EXISTS
                        )
else()
  message("      compiling without feedback file support")
  ecbuild_add_executable( TARGET  orcamodel_hofx.x
                          SOURCES orcamodelHofX.cc
                          LIBS    orcamodel
                        )

  ecbuild_add_executable( TARGET  orcamodel_hofx3D.x
                          SOURCES orcamodelHofX3D.cc
                          LIBS    orcamodel
                        )
endif()
