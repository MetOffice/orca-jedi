add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(orca-jedi)
add_subdirectory(Data)

ecbuild_add_test( TARGET test_orcajedi_coding_norms
                  TYPE SCRIPT
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../tools/cpplint.py
                  ARGS --recursive ${CMAKE_CURRENT_SOURCE_DIR}/../ )

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

ecbuild_add_test( TARGET test_orcamodel_hofx_ssh_parallel_checkerboard
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh_checkerboard.yaml
                  COMMAND orcamodel_hofx.x )

# Requires atlas-orca halos > 0 change (https://github.com/ecmwf/atlas-orca/pull/20)
ecbuild_add_test( TARGET test_orcamodel_hofx_ssh_parallel_eorca025
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh_eorca025.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_ssh_eorca025
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx_nc_ssh_eorca025.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_ssh_amm1
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx_nc_ssh_amm1.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_ssh_amm1
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx3d_nc_ssh_amm1.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_ssh_parallel_serial
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_ssh
                  ARGS testinput/hofx_nc_ssh.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_potm
                  OMP 1
                  ARGS testinput/hofx3d_nc_potm.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx3D_prof_2var
                  OMP 1
                  ARGS testinput/hofx3d_nc_prof_2vars.yaml
                  COMMAND orcamodel_hofx3D.x )
