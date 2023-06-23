add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(orca-jedi)
add_subdirectory(Data)

ecbuild_add_test( TARGET test_orcamodel_hofx_ice
                  OMP 1
                  ARGS testinput/hofx_nc_ice.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_ice
                  OMP 1
                  ARGS testinput/ostia_seaice_obs_hofx3d.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_sst
                  OMP 1
                  ARGS testinput/hofx3d_nc_sst.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_sst
                  OMP 1
                  ARGS testinput/hofx_nc_sst.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_potm
                  OMP 1
                  ARGS testinput/hofx3d_nc_potm.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_prof_2var
                  OMP 1
                  ARGS testinput/hofx3d_nc_prof_2vars.yaml
                  COMMAND orcamodel_hofx3D.x )
