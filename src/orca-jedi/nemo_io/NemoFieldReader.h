
/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <experimental/filesystem>
#include <netcdf>

#include <string>
#include <memory>
#include <vector>

#include "eckit/filesystem/PathName.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"

namespace orcamodel {

/// \brief Manager of a NEMO field file for reading the pseudo-model state
class NemoFieldReader : private util::ObjectCounter<NemoFieldReader> {
 public:
  static const std::string classname() {return "orcamodel::NemoFieldReader";}

  explicit NemoFieldReader(const eckit::PathName& filename);

  void read_datetimes();
  std::vector<atlas::PointXY> read_locs() const;
  size_t read_dim_size(const std::string& name) const;
  size_t get_nearest_datetime_index(const util::DateTime& datetime) const;
  template<typename T> T read_fillvalue(const std::string& name) const;
  template<typename T> std::vector<T> read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) const;
  template<typename T> std::vector<T> read_vertical_var(const std::string& varname,
      const size_t nlevels) const;

 private:
  NemoFieldReader() : ncFile() {}
  std::unique_ptr<netCDF::NcFile> ncFile;
  std::vector<util::DateTime> datetimes_;
  std::string time_dimvar_name_;
  std::string z_dimvar_name_;
};
}  // namespace orcamodel

//------------------------------------------------------------------------------
