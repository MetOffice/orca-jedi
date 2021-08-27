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

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

#include "atlas/field.h"

namespace orcamodel {

NemoFieldReader::NemoFieldReader( eckit::PathName& filename )
  : ncFile(nullptr)
  , datetimes_() {
  oops::Log::debug() << "orcamodel::NemoFieldReader::NemoFieldReader filename : "
                     << filename.fullName().asString() << std::endl;
  std::ostringstream err_stream;
  if( !(filename.exists()) ) {
     err_stream << "orcamodel::NemoFieldReader::NemoFieldReader filename: "
                << filename.fullName().asString() << " doesn't exist " << std::endl;
     eckit::BadValue(err_stream.str(), Here());
  }
  ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(), netCDF::NcFile::read);
  if(ncFile->isNull()) {
    err_stream << "orcamodel::NemoFieldReader::NemoFieldReader cannot open " << filename
               << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  time_dimvar_name_ = "";
  netCDF::NcDim nc_dim_time;
  netCDF::NcVar nc_var_time;
  std::vector<std::string> possible_t_dimvars({"t", "time", "time_counter"});
  for (auto & candidate_t_dimvar : possible_t_dimvars ) {
    ncFile->getCoordVar(candidate_t_dimvar, nc_dim_time, nc_var_time);
    if( (!nc_dim_time.isNull()) && (!nc_var_time.isNull()) ) {
      time_dimvar_name_ = candidate_t_dimvar;
      break;
    }
  }

  if (time_dimvar_name_ == "") {
    err_stream << "orcamodel::NemoFieldReader::read_datetimes ncVar time coordinate"
               << " is not present in NetCDF file" << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  read_datetimes( datetimes_);
}

size_t NemoFieldReader::read_dim_size( const std::string& name ) {

  auto dim = ncFile->getDim( name );
  if ( dim.isNull() ) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_dim_size Dimension '"
               << name << "' is not present in NetCDF file" << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }
  oops::Log::debug() << "orcamodel::NemoFieldReader:: group name "
                     << ncFile->getName(true) << " dim name: "  << name <<std::endl;

  return dim.getSize();
}

/// \brief retrieve the datetime for each time index in file.
void NemoFieldReader::read_datetimes(std::vector<util::DateTime>& datetimes) {

  // read time indices from file
  size_t n_times = read_dim_size(time_dimvar_name_);

  netCDF::NcVar nc_var_time = ncFile->getVar(time_dimvar_name_);
  if(nc_var_time.isNull()) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_datetimes ncVar " << time_dimvar_name_
               << " is not present in NetCDF file" << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  std::vector<int64_t> timestamps(n_times);
  nc_var_time.getVar({0}, {n_times}, timestamps.data());

  // read time units attribute from file
  netCDF::NcVarAtt nc_att_units;
  std::string units_string;
  nc_att_units = nc_var_time.getAtt("units");
  if(nc_att_units.isNull()) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_datetimes ncVar units is not present in NetCDF file";
    eckit::BadValue(err_stream.str(), Here());
  }
  nc_att_units.getValues(units_string);

  const std::string seconds_since = "seconds since ";
  std::for_each(units_string.begin(), units_string.begin()+seconds_since.size(), [](char & c){
    c = tolower(c);
  });
  if (units_string.substr(0,seconds_since.size()) != seconds_since) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_datetimes units attribute badly formatted: "
               << units_string <<std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }
  units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
  units_string.append("Z");

  auto epoch = util::DateTime(units_string.substr(seconds_since.size()));

  // construct date times
  datetimes.resize(n_times);
  for ( size_t i; i < n_times; ++i) {
    datetimes[i] = epoch + util::Duration(timestamps[i]);
  }
}

/// \brief get the time dimension index corresponding to the nearest datetime to a target datetime.
size_t NemoFieldReader::get_nearest_datetime_index(const util::DateTime& tgt_datetime) {

  int64_t time_diff = INT64_MAX;
  size_t indx = 0;
  util::Duration duration;

  for (size_t i=0; i < datetimes_.size(); ++i) {
    duration = (datetimes_[i] - tgt_datetime);
    if ( std::abs(duration.toSeconds()) < time_diff ) {
      time_diff = std::abs(duration.toSeconds());
      indx = i;
    }
  }

  return indx;
}

std::vector<atlas::PointXY> NemoFieldReader::read_locs() {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    std::string varname = "nav_lat";
    netCDF::NcVar nc_var_lat = ncFile->getVar(varname);
    if(nc_var_lat.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar " << varname
                 << " is not present in NetCDF file" << std::endl;
      eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> lats(nx*ny);

    nc_var_lat.getVar({0, 0}, {ny, nx}, lats.data());

    varname = "nav_lon";
    netCDF::NcVar nc_var_lon = ncFile->getVar(varname);
    if(nc_var_lon.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar '" << varname
                 <<"' is not present in NetCDF file.";
      eckit::BadValue(err_stream.str(), Here());
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
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_locs NetCDF exception: " <<std::endl;
    err_stream << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

std::vector<double> NemoFieldReader::read_surf_var(const std::string& varname, const size_t t_indx) {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if(nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> var_data(nx*ny);

    nc_var.getVar({t_indx, 0, 0}, {1, ny, nx}, var_data.data());

    return var_data;

  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_surf_var NetCDF exception: " <<std::endl;
    err_stream << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

void NemoFieldReader::read_surf_var(const std::string& varname, const size_t t_indx,
  atlas::array::ArrayView<double, 1>& field_view) {

  try {

    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    if (field_view.size() != nx*ny ) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var field_view dimensions "
                 << "do not match dimensions in netCDF file ";
      eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if(nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      eckit::BadValue(err_stream.str(), Here());
    }

    nc_var.getVar({t_indx, 0, 0}, {1, ny, nx}, field_view.data());


  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_surf_var NetCDF exception: " <<std::endl;
    err_stream << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}
}  // namespace orcamodel
