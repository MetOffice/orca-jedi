ecbuild_add_executable( TARGET  orcamodel_hofx.x
                        SOURCES orcamodelHofX.cc
                        LIBS    orcamodel
                      )

ecbuild_add_executable( TARGET  orcamodel_hofx3D.x
                        SOURCES orcamodelHofX3D.cc
                        LIBS    orcamodel
                      )
