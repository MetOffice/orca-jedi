add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(orca-jedi)


ecbuild_add_test( TARGET test_orcamodel_hofx3D
                  OMP 1
                  ARGS testinput/ostia_seaice_obs_hofx3d.yaml
                  COMMAND orcamodel_hofx3D.x )
