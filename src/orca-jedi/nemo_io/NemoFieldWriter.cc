/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/NemoFieldWriter.h"

#include <netcdf.h>

#include "orca-jedi/utilities/Types.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "atlas/runtime/Log.h"  // IWYU pragma: keep
#include "atlas/field.h"  // IWYU pragma: keep
#include "atlas/array.h"  // IWYU pragma: keep
#include "atlas/runtime/Exception.h"  // IWYU pragma: keep

namespace orcamodel {

namespace {

/// \brief Throw an eckit::FailedLibraryCall if a netcdf-c call failed.
void nc_check(int status, const std::string& where) {
  if (status != NC_NOERR) {
    throw eckit::FailedLibraryCall("NetCDF",
        "orcamodel::NemoFieldWriter", where + ": " + nc_strerror(status),
        Here());
  }
}

/// \brief Return true (and set varid) if a variable of the given name exists.
bool nc_var_exists(int ncid, const std::string& name, int& varid) {
  const int status = nc_inq_varid(ncid, name.c_str(), &varid);
  if (status == NC_ENOTVAR) return false;
  nc_check(status, "nc_inq_varid " + name);
  return true;
}

/// \brief Look up a dimension id by name.
int get_dimid(int ncid, const std::string& name) {
  int dimid = -1;
  nc_check(nc_inq_dimid(ncid, name.c_str(), &dimid), "nc_inq_dimid " + name);
  return dimid;
}

// Overload set dispatching a hyperslab write to the correct typed netcdf-c call.
inline int nc_put_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, const double* buf) {
  return nc_put_vara_double(ncid, varid, start, count, buf);
}
inline int nc_put_vara_dispatch(int ncid, int varid, const size_t* start,
    const size_t* count, const float* buf) {
  return nc_put_vara_float(ncid, varid, start, count, buf);
}
}  // namespace

NemoFieldWriter::NemoFieldWriter(const eckit::PathName& filename,
    const std::vector<util::DateTime>& datetimes,
    size_t nx, size_t ny,
    const std::vector<double>& depths) :
      datetimes_(datetimes)
    , depths_(depths)
    , nLevels_(depths.size())
    , nTimes_(datetimes.size())
    , nx_(nx)
    , ny_(ny)
    , dimension_variables_present_(false) {
    bool new_file;
    if (!filename.exists()) {
      nc_check(nc_create(filename.fullName().asString().c_str(),
          NC_NETCDF4 | NC_CLOBBER, &ncid_),
          "NemoFieldWriter nc_create " + filename);
      new_file = true;
    } else {
      nc_check(nc_open(filename.fullName().asString().c_str(),
          NC_WRITE, &ncid_),
          "NemoFieldWriter nc_open " + filename);
      new_file = false;
    }

    if (new_file) {
      setup_dimensions();
    }

    int varid = -1;
    dimension_variables_present_ =
        nc_var_exists(ncid_, "nav_lat", varid) &&
        nc_var_exists(ncid_, "nav_lon", varid) &&
        nc_var_exists(ncid_, "t", varid) &&
        nc_var_exists(ncid_, "z", varid);
}

NemoFieldWriter::~NemoFieldWriter() {
  if (ncid_ >= 0) {
    nc_close(ncid_);
    ncid_ = -1;
  }
}

void NemoFieldWriter::setup_dimensions() {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::setup_dimensions" << std::endl;
    int dimid = -1;
    nc_check(nc_def_dim(ncid_, "x", nx_, &dimid), "setup_dimensions def x");
    nc_check(nc_def_dim(ncid_, "y", ny_, &dimid), "setup_dimensions def y");
    nc_check(nc_def_dim(ncid_, "z", nLevels_, &dimid), "setup_dimensions def z");
    nc_check(nc_def_dim(ncid_, "t", nTimes_, &dimid), "setup_dimensions def t");
}

void NemoFieldWriter::write_dimensions(const std::vector<double>& lats,
                                       const std::vector<double>& lons) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_dimensions" << std::endl;

    if ((lons.size() != nx_*ny_) || (lats.size() != nx_*ny_)) {
      throw eckit::BadValue(
          std::string("orcamodel::NemoFieldWriter::write_dimensions")
          + " dimensions sizes do not match lat/lon buffer sizes",
        Here());
    }
    if (dimension_variables_present_) {
       oops::Log::trace() << "orcamodel::NemoFieldWriter::write_dimensions "
                          << "dimensions already present" << std::endl;
       return;
    }

    const int nxDim = get_dimid(ncid_, "x");
    const int nyDim = get_dimid(ncid_, "y");
    const int nzDim = get_dimid(ncid_, "z");
    const int ntDim = get_dimid(ncid_, "t");

    const int latlon_dims[2] = {nyDim, nxDim};
    int navLatVar = -1;
    nc_check(nc_def_var(ncid_, "nav_lat", NC_DOUBLE, 2, latlon_dims, &navLatVar),
        "write_dimensions def nav_lat");
    int navLonVar = -1;
    nc_check(nc_def_var(ncid_, "nav_lon", NC_DOUBLE, 2, latlon_dims, &navLonVar),
        "write_dimensions def nav_lon");
    {
      const size_t starts[2] = {0, 0};
      const size_t counts[2] = {ny_, nx_};
      nc_check(nc_put_vara_double(ncid_, navLonVar, starts, counts, lons.data()),
          "write_dimensions put nav_lon");
      nc_check(nc_put_vara_double(ncid_, navLatVar, starts, counts, lats.data()),
          "write_dimensions put nav_lat");
    }

    const int t_dims[1] = {ntDim};
    int timeVar = -1;
    nc_check(nc_def_var(ncid_, "t", NC_INT, 1, t_dims, &timeVar),
        "write_dimensions def t");
    const std::string seconds_since = "seconds since ";
    std::string units_string = "seconds since 1970-01-01 00:00:00";
    nc_check(nc_put_att_text(ncid_, timeVar, "units",
        units_string.size(), units_string.c_str()),
        "write_dimensions put t units");

    units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
    units_string.append("Z");
    util::DateTime epoch = util::DateTime(
        units_string.substr(seconds_since.size()));
    for (size_t iTime = 0; iTime < nTimes_; ++iTime) {
      int seconds_since_epoch = (datetimes_[iTime] - epoch).toSeconds();
      const size_t index[1] = {iTime};
      nc_check(nc_put_var1_int(ncid_, timeVar, index, &seconds_since_epoch),
          "write_dimensions put t");
    }

    {
      const int z_dims[1] = {nzDim};
      int levVar = -1;
      nc_check(nc_def_var(ncid_, "z", NC_DOUBLE, 1, z_dims, &levVar),
          "write_dimensions def z");
      const size_t starts[1] = {0};
      const size_t counts[1] = {nLevels_};
      nc_check(nc_put_vara_double(ncid_, levVar, starts, counts, depths_.data()),
          "write_dimensions put z");
    }

    dimension_variables_present_ = true;
}

template <typename T> void NemoFieldWriter::write_surf_var(const std::string& varname,
    const std::vector<T>& var_data, size_t iTime) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_surf_var" << std::endl;
    if (!dimension_variables_present_) {
      throw eckit::BadValue(
          std::string("orcamodel::NemoFieldWriter::write_surf_var")
          + " can't write '" + varname + "' as the dimensions have not yet been constructed",
        Here());
    }

    int varid = -1;
    if (!nc_var_exists(ncid_, varname, varid)) {
      const int dims[3] = {get_dimid(ncid_, "t"), get_dimid(ncid_, "y"),
                           get_dimid(ncid_, "x")};
      nc_check(nc_def_var(ncid_, varname.c_str(), NetCDFTypeMap<T>::ncType,
          3, dims, &varid), "write_surf_var def " + varname);
      const T fill = static_cast<T>(NemoFieldWriter::fillValue);
      nc_check(nc_def_var_fill(ncid_, varid, 0, &fill),
          "write_surf_var fill " + varname);
    }

    const size_t starts[3] = {iTime, 0, 0};
    const size_t counts[3] = {1, ny_, nx_};
    nc_check(nc_put_vara_dispatch(ncid_, varid, starts, counts, var_data.data()),
        "write_surf_var put " + varname);
}

template void NemoFieldWriter::write_surf_var<double>(const std::string& varname,
    const std::vector<double>& var_data, size_t iTime);
template void NemoFieldWriter::write_surf_var<float>(const std::string& varname,
    const std::vector<float>& var_data, size_t iTime);

template <typename T> void NemoFieldWriter::write_vol_var(const std::string& varname,
    const std::vector<T>& var_data, size_t iTime) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_vol_var" << std::endl;
    if (!dimension_variables_present_) {
      throw eckit::BadValue(
          std::string("orcamodel::NemoFieldWriter::write_vol_var")
          + " can't write '" + varname + "' as the dimensions have not yet been constructed",
        Here());
    }

    int varid = -1;
    if (!nc_var_exists(ncid_, varname, varid)) {
      const int dims[4] = {get_dimid(ncid_, "t"), get_dimid(ncid_, "z"),
                           get_dimid(ncid_, "y"), get_dimid(ncid_, "x")};
      nc_check(nc_def_var(ncid_, varname.c_str(), NetCDFTypeMap<T>::ncType,
          4, dims, &varid), "write_vol_var def " + varname);
      const T fill = static_cast<T>(NemoFieldWriter::fillValue);
      nc_check(nc_def_var_fill(ncid_, varid, 0, &fill),
          "write_vol_var fill " + varname);
    }

    const size_t starts[4] = {iTime, 0, 0, 0};
    const size_t counts[4] = {1, nLevels_, ny_, nx_};
    nc_check(nc_put_vara_dispatch(ncid_, varid, starts, counts, var_data.data()),
        "write_vol_var put " + varname);
}

template void NemoFieldWriter::write_vol_var<double>(const std::string& varname,
    const std::vector<double>& var_data, size_t iTime);
template void NemoFieldWriter::write_vol_var<float>(const std::string& varname,
    const std::vector<float>& var_data, size_t iTime);
}  // namespace orcamodel
