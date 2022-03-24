add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(orca-jedi)
add_subdirectory(Data)


ecbuild_add_test( TARGET test_orcamodel_hofx3D_ice
                  OMP 1
                  ARGS testinput/ostia_seaice_obs_hofx3d.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_sst
                  OMP 1
                  ARGS testinput/hofx3d_nc_sst.yaml
                  COMMAND orcamodel_hofx3D.x )
