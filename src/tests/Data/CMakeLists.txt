
list( APPEND orcajedi_test_data
hofx_sic_ostia_obs.odb
hofx_sst_obs.nc
eorca025_t_nemo.nc
orca2_t_bkg_var.nc
orca2_t_coords.nc
orca2_t_nemo.nc
orca1_t_nemo.nc
hofx_sic_ostia_obs.nc
simple_nemo.nc
create_orca2_t_bkg.nco
create_orca2_t_bkg_var.nco
hofx_sst_obs.cdl
orca2_t_coords.cdl
hofx_sic_ostia_obs.csv
simple_nemo.cdl
hofx_two_vars_obs.nc
hofx_potm_obs.nc
hofx_prof_2var_obs.nc
hofx_ssh_obs.nc
amm1_coords.nc
amm1_nemo.nc
amm1_atlas_grid_spec.yaml
amm1r_nemo.nc
amm1r_atlas_grid_spec.yaml
hofx_potm_amm1r_obs.nc
hofx_ssh_amm_obs.nc
sic_obs_ideal.nc
)

foreach(FILENAME ${orcajedi_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
