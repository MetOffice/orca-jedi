add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(orca-jedi)
add_subdirectory(Data)

# We can make use of CMake's native support of cpplint in future. If we set
# `CMAKE_<LANG>_CPPLINT` (where `<LANG>` is either `C` or `CXX`) to point to
# the `cpplint` executable, it will run `cpplint` as part of the build.
find_program(
  CPPLINT_EXECUTABLE
  NAMES cpplint
  DOC "Path to the cpplint executable"
)
if( CPPLINT_EXECUTABLE )
  ecbuild_add_test(
    TARGET "${PROJECT_NAME}_coding_norms"
    LABELS lint
    TYPE SCRIPT
    COMMAND "${CPPLINT_EXECUTABLE}"
    ARGS --quiet --recursive "${CMAKE_CURRENT_SOURCE_DIR}/.."
  )
endif()

ecbuild_add_test( TARGET orcajedi_hofx_sic
                  OMP 1
                  ARGS testinput/hofx_nc_sic.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_sic
                  OMP 1
                  ARGS testinput/hofx3d_nc_sic.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_sst
                  OMP 1
                  ARGS testinput/hofx3d_nc_sst.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_hofx_sst
                  OMP 1
                  ARGS testinput/hofx_nc_sst.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx_ssh_parallel_checkerboard
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh_checkerboard.yaml
                  COMMAND orcamodel_hofx.x )

# Requires atlas-orca halos > 0 change (https://github.com/ecmwf/atlas-orca/pull/20)
ecbuild_add_test( TARGET orcajedi_hofx_ssh_parallel_eorca025
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh_parallel_eorca025.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx_ssh_eorca025
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx_nc_ssh_eorca025.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx_ssh_amm1
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx_nc_ssh_amm1.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_ssh_amm1
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx3d_nc_ssh_amm1.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_hofx_potm_amm1r
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx_nc_potm_amm1r.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_ssh_amm1r
                  OMP 1
                  MPI 1
                  ARGS testinput/hofx3d_nc_ssh_amm1r.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_hofx_ssh_parallel_serial
                  OMP 1
                  MPI 2
                  ARGS testinput/hofx_nc_ssh_parallel.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx_ssh
                  ARGS testinput/hofx_nc_ssh.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_potm
                  OMP 1
                  ARGS testinput/hofx3d_nc_potm.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_hofx3D_prof_2var
                  OMP 1
                  ARGS testinput/hofx3d_nc_prof_2vars.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET orcajedi_3DVar_sic
                  OMP 1
                  ARGS testinput/3dvar_sic.yaml
                  COMMAND orcamodel_3DVar.x )

ecbuild_add_test( TARGET orcajedi_3DVar_sic_identity
                  OMP 1
                  ARGS testinput/3dvar_sic_identity.yaml
                  COMMAND orcamodel_3DVar.x )

ecbuild_add_test( TARGET orcajedi_Dirac
                  OMP 1
                  ARGS testinput/dirac.yaml
                  COMMAND orcamodel_ErrorCovarianceToolbox.x )
                  
ecbuild_add_test( TARGET orcajedi_bump_nicas_setup
                  OMP 1
                  ARGS testinput/bump_nicas_setup.yaml
                  COMMAND orcamodel_ErrorCovarianceToolbox.x )

