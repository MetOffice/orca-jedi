/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "NemoFieldReader.h"

#include <netcdf>
//#include <netcdf> if using Lynton Appel's netcdf-cxx4 from
//https://github.com/Unidata/netcdf-cxx4

#include <algorithm>
#include <sstream>

#include "oops/util/Logger.h"

#include "atlas/field.h"

namespace orcamodel {

NemoFieldReader::NemoFieldReader( eckit::PathName& filename ) {
  oops::Log::debug() << "orcamodel::NemoFieldReader::NemoFieldReader filename : " 
                     << filename.fullName().asString() << std::endl;
  if( !(filename.exists()) ) {
    throw("orcamodel::NemoFieldReader::NemoFieldReader filename doesn't exist ");
  }
  ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(), netCDF::NcFile::read);
  if(ncFile->isNull()) {
    throw("orcamodel::NemoFieldReader::NemoFieldReader " + filename + " not found");
  }
}

size_t NemoFieldReader::read_dim_size( const std::string& name ) {

  auto dim = ncFile->getDim( name );
  if ( dim.isNull() ) {
    throw ( "orcamodel::NemoFieldReader::read_locs Dimension '" + name + "' is not present in NetCDF file" );
  }
  oops::Log::debug() << "orcamodel::NemoFieldReader:: group name " 
                     << ncFile->getName(true) << " dim name: "  << name <<std::endl;

  return dim.getSize();
}

std::vector<atlas::PointXY> NemoFieldReader::read_locs() {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    std::string varname = "nav_lat";
    netCDF::NcVar nc_var_lat = ncFile->getVar(varname);
    if(nc_var_lat.isNull()) {
      throw("orcamodel::NemoFieldReader::read_locs ncVar " + varname + " is not present in NetCDF file");
    }

    std::vector<double> lats(nx*ny);

    nc_var_lat.getVar({0, 0}, {ny, nx}, lats.data());

    varname = "nav_lon";
    netCDF::NcVar nc_var_lon = ncFile->getVar(varname);
    if(nc_var_lon.isNull()) {
      throw("orcamodel::NemoFieldReader::read_locs ncVar " + varname + " is not present in NetCDF file");
    }

    std::vector<double> lons(nx*ny);

    nc_var_lon.getVar({0, 0}, {ny, nx}, lons.data());

    std::vector<atlas::PointXY> locations(nx*ny);
    
    for (size_t i=0; i<nx*ny; i++){
      locations[i] = atlas::PointXY(lons[i], lats[i]);
    }

    return locations;

  } catch(netCDF::exceptions::NcException& e)
  {
    e.what();
    oops::Log::debug() << "orcamodel::NemoFieldReader::read_locs ERROR: " <<std::endl;
    throw(e);
  }
}

std::vector<double> NemoFieldReader::read_surf_var(std::string varname) {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if(nc_var.isNull()) {
      throw("orcamodel::NemoFieldReader::read_surf_var ncVar " 
            + varname + " is not present in NetCDF file");
    }

    std::vector<double> var_data(nx*ny);

    nc_var.getVar({0, 0, 0}, {1, ny, nx}, var_data.data());

    return var_data;

  } catch(netCDF::exceptions::NcException& e)
  {
    e.what();
    oops::Log::debug() << "orcamodel::NemoFieldReader::read_surf_var ERROR: " <<std::endl;
    throw(e);
  }
}

void NemoFieldReader::read_surf_var(std::string varname, atlas::array::ArrayView<double, 1>& field_view) {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    if (field_view.size() != nx*ny ) {
      throw("orcamodel::NemoFieldReader::read_surf_var field_view dimensions "
            "do not match dimensions in netCDF file ");
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if(nc_var.isNull()) {
      throw("orcamodel::NemoFieldReader::read_surf_var ncVar " 
            + varname + " is not present in NetCDF file");
    }

    nc_var.getVar({0, 0, 0}, {1, ny, nx}, field_view.data());


  } catch(netCDF::exceptions::NcException& e)
  {
    e.what();
    oops::Log::debug() <<  "orcamodel::NemoFieldReader::read_surf_var ERROR: " <<std::endl;
    throw(e);
  }
}
}  // namespace orcamodel
