/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include <netcdf.h>

#include <algorithm>
#include <sstream>
#include <limits>
#include <utility>
#include <vector>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

namespace orcamodel {

namespace {

/// \brief Throw an eckit::ReadError if a netcdf-c call returned an error.
void nc_check(int status, const std::string& where) {
  if (status != NC_NOERR) {
    throw eckit::ReadError(std::string("orcamodel::NemoFieldReader ") + where
        + ": " + nc_strerror(status), Here());
  }
}

/// \brief Return true (and set varid) if a variable of the given name exists.
bool nc_var_exists(int ncid, const std::string& name, int& varid) {
  const int status = nc_inq_varid(ncid, name.c_str(), &varid);
  if (status == NC_ENOTVAR) return false;
  nc_check(status, "nc_inq_varid " + name);
  return true;
}

/// \brief Return true (and set dimid) if a dimension of the given name exists.
bool nc_dim_exists(int ncid, const std::string& name, int& dimid) {
  const int status = nc_inq_dimid(ncid, name.c_str(), &dimid);
  if (status == NC_EBADDIM) return false;
  nc_check(status, "nc_inq_dimid " + name);
  return true;
}

/// \brief Look up a variable id, throwing eckit::BadValue if it is absent.
int get_varid_or_throw(int ncid, const std::string& name,
    const std::string& where) {
  int varid = -1;
  if (!nc_var_exists(ncid, name, varid)) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::" << where << " ncVar '"
               << name << "' is not present in NetCDF file";
    throw eckit::BadValue(err_stream.str(), Here());
  }
  return varid;
}

/// \brief Search the netCDF file for a dimension matching a name from a list of
///        possible names.
/// \param ncid The netCDF file id.
/// \param check_dim_for_dimvar Set to true to ensure the dimension has
///        a corresponding dimension variable.
/// \param possible_names A vector of all possible names for the variable.
/// \return The dimension name.
std::string find_nc_var_name(int ncid, const bool check_dim_for_dimvar,
    const std::vector<std::string>& possible_names) {
  for (const auto& candidate : possible_names) {
    int dimid = -1;
    int varid = -1;
    const bool has_dim = nc_dim_exists(ncid, candidate, dimid);
    const bool has_var = nc_var_exists(ncid, candidate, varid);
    if (has_dim && (!check_dim_for_dimvar || has_var)) {
      return candidate;
    }
  }

  std::ostringstream err_stream;
  err_stream << "orcamodel::find_nc_var_name coordinate matching {";
  for (const auto& n : possible_names) err_stream << "'" << n << "' ";
  err_stream << "} is not present in NetCDF file" << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

// Overload set dispatching a hyperslab read to the correct typed netcdf-c call.
inline int nc_get_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, double* buf) {
  return nc_get_vara_double(ncid, varid, start, count, buf);
}
inline int nc_get_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, float* buf) {
  return nc_get_vara_float(ncid, varid, start, count, buf);
}
inline int nc_get_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, int* buf) {
  return nc_get_vara_int(ncid, varid, start, count, buf);
}
inline int nc_get_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, long long* buf) {  // NOLINT(runtime/int)
  return nc_get_vara_longlong(ncid, varid, start, count, buf);
}

/// \brief Read a hyperslab of a variable stored as InType, casting to OutType.
template<typename InType, typename OutType>
void fill_from_ncvar_as_type(int ncid, int varid,
    const std::vector<size_t>& starts, const std::vector<size_t>& counts,
    std::vector<OutType>& result) {
  std::vector<InType> buffer(result.size());
  nc_check(nc_get_vara_dispatch(ncid, varid, starts.data(), counts.data(),
      buffer.data()), "fill_from_ncvar_as_type nc_get_vara");
  std::transform(buffer.begin(), buffer.end(), result.begin(),
     [](const InType in) -> OutType { return static_cast<OutType>(in); });
}

template<>
void fill_from_ncvar_as_type<double, double>(int ncid, int varid,
    const std::vector<size_t>& starts, const std::vector<size_t>& counts,
    std::vector<double>& result) {
  nc_check(nc_get_vara_double(ncid, varid, starts.data(), counts.data(),
      result.data()), "fill_from_ncvar_as_type nc_get_vara_double");
}
template<>
void fill_from_ncvar_as_type<float, float>(int ncid, int varid,
    const std::vector<size_t>& starts, const std::vector<size_t>& counts,
    std::vector<float>& result) {
  nc_check(nc_get_vara_float(ncid, varid, starts.data(), counts.data(),
      result.data()), "fill_from_ncvar_as_type nc_get_vara_float");
}

// Overload set dispatching a scalar attribute read to the correct typed call.
inline int nc_get_att_dispatch(int ncid, int varid, const char* name,
    double* v) {
  return nc_get_att_double(ncid, varid, name, v);
}
inline int nc_get_att_dispatch(int ncid, int varid, const char* name,
    float* v) {
  return nc_get_att_float(ncid, varid, name, v);
}
inline int nc_get_att_dispatch(int ncid, int varid, const char* name,
    int* v) {
  return nc_get_att_int(ncid, varid, name, v);
}
}  // namespace

NemoFieldReader::NemoFieldReader(const eckit::PathName& filename)
  : ncid_(-1), datetimes_() {
  oops::Log::debug() << "orcamodel::NemoFieldReader::NemoFieldReader filename: "
                     << filename.fullName().asString() << std::endl;
  if (!filename.exists()) {
     std::ostringstream err_stream;
     err_stream << "orcamodel::NemoFieldReader::NemoFieldReader filename: "
                << filename.fullName().asString() << " doesn't exist "
                << std::endl;
     throw eckit::BadValue(err_stream.str(), Here());
  }
  // netcdf-c's nc_open reads netCDF-3 (classic) and netCDF-4 (HDF5)
  // transparently, so the serial read path supports any on-disk format.
  const int status = nc_open(filename.fullName().asString().c_str(),
                             NC_NOWRITE, &ncid_);
  if (status != NC_NOERR) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::NemoFieldReader cannot open "
               << filename << ": " << nc_strerror(status) << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  time_dimvar_name_ = find_nc_var_name(ncid_, true,
                                       {"t", "time", "time_counter"});
  z_dimvar_name_ = find_nc_var_name(ncid_, false, {"z", "deptht"});

  read_datetimes();
}

NemoFieldReader::~NemoFieldReader() {
  if (ncid_ >= 0) {
    nc_close(ncid_);
  }
}

NemoFieldReader::NemoFieldReader(NemoFieldReader&& other) noexcept
  : ncid_(other.ncid_),
    datetimes_(std::move(other.datetimes_)),
    time_dimvar_name_(std::move(other.time_dimvar_name_)),
    z_dimvar_name_(std::move(other.z_dimvar_name_)) {
  other.ncid_ = -1;
}

NemoFieldReader& NemoFieldReader::operator=(NemoFieldReader&& other) noexcept {
  if (this != &other) {
    if (ncid_ >= 0) nc_close(ncid_);
    ncid_ = other.ncid_;
    datetimes_ = std::move(other.datetimes_);
    time_dimvar_name_ = std::move(other.time_dimvar_name_);
    z_dimvar_name_ = std::move(other.z_dimvar_name_);
    other.ncid_ = -1;
  }
  return *this;
}

/// \brief Read the dimension size for a given netCDF dimension specified by name
/// \param name The name of the netCDF dimension
/// \return size The size
size_t NemoFieldReader::read_dim_size(const std::string& name) const {
  int dimid = -1;
  if (!nc_dim_exists(ncid_, name, dimid)) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_dim_size Dimension '"
               << name << "' is not present in NetCDF file" << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  size_t len = 0;
  nc_check(nc_inq_dimlen(ncid_, dimid, &len), "read_dim_size " + name);
  return len;
}

/// \brief Update the datetimes_ in the object to contain times for each time index in file.
void NemoFieldReader::read_datetimes() {
  // read time indices from file
  size_t n_times = read_dim_size(time_dimvar_name_);
  const int varid = get_varid_or_throw(ncid_, time_dimvar_name_,
                                       "read_datetimes");

  if (n_times < 1) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_datetimes n_times < 1 "
               << n_times << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  std::vector<long long> timestamps(n_times);  // NOLINT(runtime/int)
  const size_t start = 0;
  nc_check(nc_get_vara_longlong(ncid_, varid, &start, &n_times,
      timestamps.data()), "read_datetimes get time values");

  // read time units attribute from file
  size_t units_len = 0;
  nc_check(nc_inq_attlen(ncid_, varid, "units", &units_len),
      "read_datetimes units attribute length");
  std::string units_string(units_len, '\0');
  nc_check(nc_get_att_text(ncid_, varid, "units", &units_string[0]),
      "read_datetimes get units attribute");

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
}

/// \brief Read the _FillValue for a variable, defaulting to the minimum value
/// for the datatype.
/// \param name Name of the netCDF variable containing the _FillValue attribute to retrieve.
/// \return The fill value for this netCDF variable
template<typename T> T NemoFieldReader::read_fillvalue(const std::string& name) const {
  const int varid = get_varid_or_throw(ncid_, name, "read_fillvalue");

  T fillvalue = std::numeric_limits<T>::lowest();

  const int status = nc_inq_att(ncid_, varid, "_FillValue", nullptr, nullptr);
  if (status == NC_ENOTATT) {
    oops::Log::trace() << "orcamodel::NemoFieldReader::read_fillvalue fillvalue"
                       << " not found " << std::endl;
  } else {
    nc_check(status, "read_fillvalue inq _FillValue " + name);
    oops::Log::trace() << "orcamodel::NemoFieldReader::read_fillvalue found: "
                       << "_FillValue" << std::endl;
    nc_check(nc_get_att_dispatch(ncid_, varid, "_FillValue", &fillvalue),
        "read_fillvalue get _FillValue " + name);
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
  size_t nx = read_dim_size("x");
  size_t ny = read_dim_size("y");

  const int lat_id = get_varid_or_throw(ncid_, "nav_lat", "read_locs");
  nc_type lat_type;
  nc_check(nc_inq_vartype(ncid_, lat_id, &lat_type), "read_locs nav_lat type");
  std::vector<double> lats(nx*ny);
  if (lat_type == NC_DOUBLE) {
    fill_from_ncvar_as_type<double, double>(ncid_, lat_id, {0, 0}, {ny, nx}, lats);
  } else if (lat_type == NC_FLOAT) {
    fill_from_ncvar_as_type<float, double>(ncid_, lat_id, {0, 0}, {ny, nx}, lats);
  } else {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_locs ncVar 'nav_lat'"
               << " reading type " << lat_type << " not supported.";
    throw eckit::BadValue(err_stream.str(), Here());
  }

  const int lon_id = get_varid_or_throw(ncid_, "nav_lon", "read_locs");
  nc_type lon_type;
  nc_check(nc_inq_vartype(ncid_, lon_id, &lon_type), "read_locs nav_lon type");
  std::vector<double> lons(nx*ny);
  if (lon_type == NC_DOUBLE) {
    fill_from_ncvar_as_type<double, double>(ncid_, lon_id, {0, 0}, {ny, nx}, lons);
  } else if (lon_type == NC_FLOAT) {
    fill_from_ncvar_as_type<float, double>(ncid_, lon_id, {0, 0}, {ny, nx}, lons);
  } else {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_locs ncVar 'nav_lon'"
               << " reading type " << lon_type << " not supported.";
    throw eckit::BadValue(err_stream.str(), Here());
  }

  std::vector<atlas::PointXY> locations(nx*ny);

  for (size_t i=0; i < nx*ny; i++) {
    locations[i] = atlas::PointXY(lons[i], lats[i]);
  }

  return locations;
}

/// \brief Read all data for a given variable at a level and time.
/// \param varname The name of the variable.
/// \param t_indx The time index of the slice.
/// \param z_indx The vertical index of the slice.
/// \return a vector of the variable data.
template<typename T> std::vector<T> NemoFieldReader::read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) const {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_var_slice("
                     << varname << ", " << t_indx << ", " << z_indx
                     << ")" << std::endl;
  size_t nx = read_dim_size("x");
  size_t ny = read_dim_size("y");

  const int varid = get_varid_or_throw(ncid_, varname, "read_var_slice");

  int n_dims = 0;
  nc_check(nc_inq_varndims(ncid_, varid, &n_dims), "read_var_slice ndims " + varname);
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

  std::vector<T> var_data(nx*ny);

  nc_type xtype;
  nc_check(nc_inq_vartype(ncid_, varid, &xtype), "read_var_slice type " + varname);
  if (xtype == NC_DOUBLE) {
      fill_from_ncvar_as_type<double, T>(ncid_, varid, starts, counts, var_data);
  } else if (xtype == NC_FLOAT) {
      fill_from_ncvar_as_type<float, T>(ncid_, varid, starts, counts, var_data);
  } else if (xtype == NC_INT) {
      fill_from_ncvar_as_type<int, T>(ncid_, varid, starts, counts, var_data);
  } else if (xtype == NC_INT64) {
      // NOLINTNEXTLINE(runtime/int)
      fill_from_ncvar_as_type<long long, T>(ncid_, varid, starts, counts, var_data);
  } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_var_slice ncVar '"
                 << varname << "' reading type "
                 << xtype << " not supported.";
      throw eckit::BadValue(err_stream.str(), Here());
  }

  return var_data;
}
template std::vector<double> NemoFieldReader::read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) const;
template std::vector<float> NemoFieldReader::read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) const;

/// \brief Read a 1D variable containing level depth information.
/// \param varname The name of the variable.
/// \param n_levels The number of levels to read (beginning from the surface).
/// \return a vector of the depth values.
template<typename T> std::vector<T> NemoFieldReader::read_vertical_var(
    const std::string& varname,
    const size_t nlevels) const {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_vertical_var"
                     << std::endl;
  size_t nz = read_dim_size(z_dimvar_name_);

  if (nlevels > nz) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_vertical_var num levels "
               << "is larger than NetCDF file z dimension ";
    throw eckit::BadValue(err_stream.str(), Here());
  }

  const int varid = get_varid_or_throw(ncid_, varname, "read_vertical_var");

  int n_dims = 0;
  nc_check(nc_inq_varndims(ncid_, varid, &n_dims),
      "read_vertical_var ndims " + varname);
  std::string first_dim_name;
  if (n_dims == 1) {
    int dimids[1] = {-1};
    nc_check(nc_inq_vardimid(ncid_, varid, dimids),
        "read_vertical_var vardimid " + varname);
    char dim_name[NC_MAX_NAME + 1] = {0};
    nc_check(nc_inq_dimname(ncid_, dimids[0], dim_name),
        "read_vertical_var dimname " + varname);
    first_dim_name = dim_name;
  }
  if (n_dims != 1 || first_dim_name != z_dimvar_name_) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
               << varname << "' has " << n_dims << " dimensions.";
    throw eckit::BadValue(err_stream.str(), Here());
  }

  std::vector<T> buffer(nlevels);
  std::vector<size_t> starts = {0, };
  std::vector<size_t> counts = {nlevels, };
  nc_type xtype;
  nc_check(nc_inq_vartype(ncid_, varid, &xtype),
      "read_vertical_var type " + varname);
  if (xtype == NC_DOUBLE) {
      fill_from_ncvar_as_type<double, T>(ncid_, varid, starts, counts, buffer);
  } else if (xtype == NC_FLOAT) {
      fill_from_ncvar_as_type<float, T>(ncid_, varid, starts, counts, buffer);
  } else if (xtype == NC_INT) {
      fill_from_ncvar_as_type<int, T>(ncid_, varid, starts, counts, buffer);
  } else if (xtype == NC_INT64) {
      // NOLINTNEXTLINE(runtime/int)
      fill_from_ncvar_as_type<long long, T>(ncid_, varid, starts, counts, buffer);
  } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' reading type "
                 << xtype << " not supported.";
      throw eckit::BadValue(err_stream.str(), Here());
  }

  return buffer;
}
template std::vector<double> NemoFieldReader::read_vertical_var(
    const std::string& varname,
    const size_t nlevels) const;
template std::vector<float> NemoFieldReader::read_vertical_var(
    const std::string& varname,
    const size_t nlevels) const;

}  // namespace orcamodel
