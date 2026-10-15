/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <netcdf.h>

#include <string>

#include "eckit/exception/Exceptions.h"

namespace orcamodel {

/// \brief Throw an eckit::FailedLibraryCall if a netcdf-c call failed.
inline void nc_check(int status, const std::string& where) {
  if (status != NC_NOERR) {
    throw eckit::FailedLibraryCall("NetCDF", where, nc_strerror(status), Here());
  }
}

/// \brief Return true (and set varid) if a variable of the given name exists.
inline bool nc_var_exists(int ncid, const std::string& name, int& varid) {
  const int status = nc_inq_varid(ncid, name.c_str(), &varid);
  if (status == NC_ENOTVAR) return false;
  nc_check(status, "nc_inq_varid " + name);
  return true;
}

/// \brief Return true (and set dimid) if a dimension of the given name exists.
inline bool nc_dim_exists(int ncid, const std::string& name, int& dimid) {
  const int status = nc_inq_dimid(ncid, name.c_str(), &dimid);
  if (status == NC_EBADDIM) return false;
  nc_check(status, "nc_inq_dimid " + name);
  return true;
}

}  // namespace orcamodel
