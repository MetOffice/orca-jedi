
list( APPEND orcajedi_test_output
  test_hofx_nc_sic.ref
  test_hofx_nc_sst.ref
  test_hofx_nc_ssh.ref
  test_hofx_nc_ssh_eorca025.ref
  test_hofx_nc_ssh_amm1.ref
  test_hofx_nc_potm_amm1r.ref
  test_hofx3d_nc_sic.ref
  test_hofx3d_nc_sst.ref
  test_hofx3d_nc_potm.ref
  test_hofx3d_nc_prof_2vars.ref
  test_hofx3d_nc_ssh_amm1.ref
  test_hofx3d_nc_ssh_amm1r.ref
  test_3dvar_sic.ref
  test_3dvar_sic_nv.ref
  test_3dvar_sic_identity.ref
  test_dirac.ref
  test_dirac_nv.ref
  test_bump_nicas.ref
)

foreach(FILENAME ${orcajedi_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
