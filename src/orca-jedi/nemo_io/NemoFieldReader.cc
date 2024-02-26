/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include <netcdf>
// Using Lynton Appel's netcdf-cxx4 from
// https://github.com/Unidata/netcdf-cxx4
//
#include <algorithm>
#include <sstream>
#include <limits>
#include <map>

#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/system/ResourceUsage.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

#include "atlas-orca/grid/OrcaGrid.h"
#include "orca-jedi/nemo_io/OrcaIndex.h"

namespace orcamodel {

namespace {

/// \brief Search the netCDF file for a dimension matching a name from a list of possible names.
/// \param ncFile The netCDF file.
/// \param check_dim_for_dimvar Set to true to ensure the dimension has
///        a corresponding dimension variable.
/// \param possible_names A vector of all possible names for the variable.
/// \return The dimension name.
std::string find_nc_var_name(const netCDF::NcFile& ncFile,
    const bool check_dim_for_dimvar,
    const std::vector<std::string>& possible_names) {
  std::string dimvar_name("");
  netCDF::NcDim nc_dim;
  netCDF::NcVar nc_var;
  for (const auto & c_candidate : possible_names) {
    // Copy because ncFile methods only have non-const interface
    std::string candidate(c_candidate);
    nc_var = ncFile.getVar(candidate);
    nc_dim = ncFile.getDim(candidate);
    // check for coordinate dimension and optionally corresponding dimension
    // variable
    if (!nc_dim.isNull() && (!check_dim_for_dimvar || !nc_var.isNull())) {
      dimvar_name = candidate;
      return dimvar_name;
    }
  }

  if (dimvar_name == "") {
    std::ostringstream err_stream;
    err_stream << "orcamodel::find_nc_var_name coordinate matching {";
    for (const auto & n : possible_names) err_stream << "'" << n << "' ";
    err_stream << "} is not present in NetCDF file" << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  return dimvar_name;
}
template<typename T> void fill_from_ncvar_as_double(
    const netCDF::NcVar nc_var,
    const std::vector<size_t>& sizes, const std::vector<size_t>& counts,
    std::vector<double>& result) {
  std::vector<T> buffer(result.size());
  nc_var.getVar(sizes, counts, buffer.data());
  std::transform(buffer.begin(), buffer.end(), result.begin(),
     [](const T in) -> double { return static_cast<double>(in); });
}
template void fill_from_ncvar_as_double<float>(
    const netCDF::NcVar nc_var,
    const std::vector<size_t>& sizes, const std::vector<size_t>& counts,
    std::vector<double>& result);
template void fill_from_ncvar_as_double<int>(
    const netCDF::NcVar nc_var,
    const std::vector<size_t>& sizes, const std::vector<size_t>& counts,
    std::vector<double>& result);
template void fill_from_ncvar_as_double<int64_t>(
    const netCDF::NcVar nc_var,
    const std::vector<size_t>& sizes, const std::vector<size_t>& counts,
    std::vector<double>& result);

template<> void fill_from_ncvar_as_double<double>(
    const netCDF::NcVar nc_var,
    const std::vector<size_t>& sizes, const std::vector<size_t>& counts,
    std::vector<double>& result) {
  nc_var.getVar(sizes, counts, result.data());
}
}  // namespace

NemoFieldReader::NemoFieldReader(const eckit::PathName& filename)
  : ncFile(nullptr), datetimes_() {
  oops::Log::debug() << "orcamodel::NemoFieldReader::NemoFieldReader filename: "
                     << filename.fullName().asString() << std::endl;
  std::ostringstream err_stream;
  if (!filename.exists()) {
     err_stream << "orcamodel::NemoFieldReader::NemoFieldReader filename: "
                << filename.fullName().asString() << " doesn't exist "
                << std::endl;
     throw eckit::BadValue(err_stream.str(), Here());
  }
  try {
    ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
        netCDF::NcFile::read);
  } catch(netCDF::exceptions::NcException& e) {
    err_stream << "orcamodel::NemoFieldReader::NemoFieldReader cannot open "
               << filename << std::endl;
    err_stream << "NetCDF Exception: " << e.what() << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  time_dimvar_name_ = find_nc_var_name(*ncFile, true,
                                       {"t", "time", "time_counter"});
  z_dimvar_name_ = find_nc_var_name(*ncFile, false, {"z", "deptht"});

  read_datetimes();
}

/// \brief Read the dimension size for a given netCDF dimension specified by name
/// \param name The name of the netCDF dimension
/// \return size The size
size_t NemoFieldReader::read_dim_size(const std::string& name) const {
  auto dim = ncFile->getDim(name);
  if (dim.isNull()) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_dim_size Dimension '"
               << name << "' is not present in NetCDF file" << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  return dim.getSize();
}

/// \brief Update the datetimes_ in the object to contain times for each time index in file.
void NemoFieldReader::read_datetimes() {
  // read time indices from file
  size_t n_times = read_dim_size(time_dimvar_name_);

  try {
    netCDF::NcVar nc_var_time = ncFile->getVar(time_dimvar_name_);
    if (nc_var_time.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_datetimes ncVar "
                 << time_dimvar_name_
                 << " is not present in NetCDF file" << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    if (n_times < 1) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_datetimes n_times < 1 "
                 << n_times << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<int64_t> timestamps(n_times);
    nc_var_time.getVar({0}, {n_times}, timestamps.data());

    // read time units attribute from file
    netCDF::NcVarAtt nc_att_units;
    std::string units_string;
    nc_att_units = nc_var_time.getAtt("units");
    if (nc_att_units.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_datetimes ncVar units is "
                 << "not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }
    nc_att_units.getValues(units_string);

    const std::string seconds_since = "seconds since ";
    std::for_each(units_string.begin(),
        units_string.begin()+seconds_since.size(),
        [](char & c){ c = tolower(c); });
    if (units_string.substr(0, seconds_since.size()) != seconds_since) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_datetimes units attribute"
                 << " badly formatted: " << units_string << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }
    units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
    units_string.append("Z");

    auto epoch = util::DateTime(units_string.substr(seconds_since.size()));

    // construct date times
    datetimes_.resize(n_times);
    for (size_t i=0; i < n_times; ++i) {
      datetimes_[i] = epoch + util::Duration(timestamps[i]);
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_datetimes NetCDF exception "
               << time_dimvar_name_ << ": " << std::endl;
    err_stream << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

/// \brief Read the _FillValue for a variable, defaulting to the minimum value
/// for the datatype.
/// \param name Name of the netCDF variable containing the _FillValue attribute to retrieve.
/// \return The fill value for this netCDF variable
template<typename T> T NemoFieldReader::read_fillvalue(const std::string& name) const {
  netCDF::NcVar nc_var = ncFile->getVar(name);
  if (nc_var.isNull()) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_fillvalue ncVar "
               << time_dimvar_name_ << " is not present in NetCDF file"
               << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  T fillvalue = std::numeric_limits<T>::lowest();

  std::map<std::string, netCDF::NcVarAtt> attributeList = nc_var.getAtts();
  auto myIter = attributeList.find("_FillValue");
  if (myIter == attributeList.end()) {
    oops::Log::trace() << "orcamodel::NemoFieldReader::read_fillvalue fillvalue"
                       << " not found " << std::endl;
  } else {
    oops::Log::trace() << "orcamodel::NemoFieldReader::read_fillvalue found: "
                       << myIter->first << std::endl;
    netCDF::NcVarAtt nc_att_fill = myIter->second;
    nc_att_fill.getValues(&fillvalue);
  }

  oops::Log::trace() << "orcamodel::NemoFieldReader::read_fillvalue fillvalue: "
                     << fillvalue << std::endl;

  return fillvalue;
}
template int NemoFieldReader::read_fillvalue<int>(const std::string& name) const;
template float NemoFieldReader::read_fillvalue<float>(const std::string& name) const;
template double NemoFieldReader::read_fillvalue<double>(
    const std::string& name) const;

/// \brief get the time dimension index corresponding to the nearest datetime
/// to a target datetime.
/// \param tgt_datetime Search for the index of the time slice in the file nearest this datetime
/// \return The index of the nearest time slice in the file
size_t NemoFieldReader::get_nearest_datetime_index(
    const util::DateTime& tgt_datetime) const {
  int64_t time_diff = INT64_MAX;
  size_t indx = 0;
  util::Duration duration;

  for (size_t i=0; i < datetimes_.size(); ++i) {
    duration = datetimes_[i] - tgt_datetime;
    if ( std::abs(duration.toSeconds()) < time_diff ) {
      time_diff = std::abs(duration.toSeconds());
      indx = i;
    }
  }

  return indx;
}

/// \brief Read the latitude longitude locations of all points in a NEMO field file
/// \return A vector of XY points of all nodes in the field.
std::vector<atlas::PointXY> NemoFieldReader::read_locs() const {
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    std::string varname = "nav_lat";
    netCDF::NcVar nc_var_lat = ncFile->getVar(varname);
    if (nc_var_lat.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar " << varname
                 << " is not present in NetCDF file" << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> lats(nx*ny);
    std::string lat_type_name = nc_var_lat.getType().getName();
    if (lat_type_name == "double") {
      fill_from_ncvar_as_double<double>(nc_var_lat, {0, 0}, {ny, nx}, lats);
    } else if (lat_type_name == "float") {
      fill_from_ncvar_as_double<float>(nc_var_lat, {0, 0}, {ny, nx}, lats);
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar '"
                 << varname << "' reading type "
                 << lat_type_name << " not supported.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    varname = "nav_lon";
    netCDF::NcVar nc_var_lon = ncFile->getVar(varname);
    if (nc_var_lon.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar '" << varname
                 <<"' is not present in NetCDF file.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> lons(nx*ny);
    std::string lon_type_name = nc_var_lat.getType().getName();
    if (lon_type_name == "double") {
      fill_from_ncvar_as_double<double>(nc_var_lon, {0, 0}, {ny, nx}, lons);
    } else if (lon_type_name == "float") {
      fill_from_ncvar_as_double<float>(nc_var_lon, {0, 0}, {ny, nx}, lons);
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar '"
                 << varname << "' reading type "
                 << nc_var_lon.getType().getName() << " not supported.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<atlas::PointXY> locations(nx*ny);

    for (size_t i=0; i < nx*ny; i++) {
      locations[i] = atlas::PointXY(lons[i], lats[i]);
    }

    return locations;
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_locs NetCDF exception: "
               << std::endl;
    err_stream << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

/// \brief Read all data for a given variable at a level and time.
/// \param varname The name of the variable.
/// \param t_indx The time index of the slice.
/// \param z_indx The vertical index of the slice.
/// \return a vector of the variable data.
std::vector<double> NemoFieldReader::read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) const {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_var_slice("
                     << varname << ", " << t_indx << ", " << z_indx
                     << ")" << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_var_slice ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    size_t n_dims = nc_var.getDimCount();
  oops::Log::trace() << "orcamodel:: read var slice " << varname << " "
      << z_indx << " before memory: "
      << static_cast<double>(eckit::system::ResourceUsage().maxResidentSetSize()) / 1.0e+6
      << " Mb" << std::endl;
    std::vector<size_t> starts;
    std::vector<size_t> counts;
    if (n_dims == 4) {
      starts = {t_indx, z_indx, 0, 0};
      counts = {1, 1, ny, nx};
    } else if (n_dims == 3) {
      starts = {t_indx, 0, 0};
      counts = {1, ny, nx};
    } else if (n_dims == 2) {
      starts = {0, 0};
      counts = {ny, nx};
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_var_slice ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> var_data(nx*ny);

    std::string nc_type_name = nc_var.getType().getName();
    if (nc_type_name == "double") {
        fill_from_ncvar_as_double<double>(nc_var, starts, counts, var_data);
    } else if (nc_type_name == "float") {
        fill_from_ncvar_as_double<float>(nc_var, starts, counts, var_data);
    } else if (nc_type_name == "int") {
        fill_from_ncvar_as_double<int>(nc_var, starts, counts, var_data);
    } else if (nc_type_name == "int64") {
        fill_from_ncvar_as_double<int64_t>(nc_var, starts, counts, var_data);
    } else {
        std::ostringstream err_stream;
        err_stream << "orcamodel::NemoFieldReader::read_var_slice ncVar '"
                   << varname << "' reading type "
                   << nc_type_name << " not supported.";
        throw eckit::BadValue(err_stream.str(), Here());
    }

    oops::Log::trace() << "orcamodel:: read var slice " << varname << " "
        << z_indx << " after vec mem " << var_data.capacity()*sizeof(double) / 1.0e+6 << " memory: "
        << static_cast<double>(eckit::system::ResourceUsage().maxResidentSetSize()) / 1.0e+6
        << " Mb" << std::endl;

    return var_data;
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_var_slice varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

/// \brief Read a 1D variable containing level depth information.
/// \param varname The name of the variable.
/// \param n_levels The number of levels to read (beginning from the surface).
/// \return a vector of the depth values.
std::vector<double> NemoFieldReader::read_vertical_var(
    const std::string& varname,
    const size_t nlevels) const {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_vertical_var"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");
    size_t nz = read_dim_size(z_dimvar_name_);

    if (nlevels > nz) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var num levels "
                 << "is larger than NetCDF file z dimension ";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    size_t n_dims = nc_var.getDimCount();
    std::string first_dim_name = nc_var.getDim(0).getName();
    if (n_dims != 1 || first_dim_name != z_dimvar_name_) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> buffer(nlevels);
    std::vector<size_t> starts = {0, };
    std::vector<size_t> counts = {nlevels, };
    std::string nc_type_name = nc_var.getType().getName();
    if (nc_type_name == "double") {
        fill_from_ncvar_as_double<double>(nc_var, starts, counts, buffer);
    } else if (nc_type_name == "float") {
        fill_from_ncvar_as_double<float>(nc_var, starts, counts, buffer);
    } else if (nc_type_name == "int") {
        fill_from_ncvar_as_double<int>(nc_var, starts, counts, buffer);
    } else if (nc_type_name == "int64") {
        fill_from_ncvar_as_double<int64_t>(nc_var, starts, counts, buffer);
    } else {
        std::ostringstream err_stream;
        err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                   << varname << "' reading type "
                   << nc_var.getType().getName() << " not supported.";
        throw eckit::BadValue(err_stream.str(), Here());
    }

    return buffer;
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_vertical_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

}  // namespace orcamodel
