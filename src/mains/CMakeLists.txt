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

  ecbuild_add_executable( TARGET  orcamodel_3DVar.x
                          SOURCES orcamodel3DVar.cc
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

  ecbuild_add_executable( TARGET  orcamodel_3DVar.x
                          SOURCES orcamodel3DVar.cc
                          LIBS    orcamodel
                        )
endif()

ecbuild_add_executable( TARGET  orcamodel_ErrorCovarianceToolbox.x
                        SOURCES orcamodelErrorCovarianceToolbox.cc
                        LIBS    orcamodel
                      )

oops_output_json_schema( "orcamodel_hofx.x" )
oops_output_json_schema( "orcamodel_hofx3D.x" )
#oops_output_json_schema( "orcamodel_3DVar.x" )
#oops_output_json_schema( "orcamodel_ErrorCovarianceToolbox.x" )
