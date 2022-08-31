
/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
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
#include "atlas/field.h"
#include "atlas/mesh.h"
#include "atlas/array/ArrayView.h"

namespace orcamodel {

/// \brief Manager of a NEMO field file for reading the pseudo-model state
class NemoFieldReader : private util::ObjectCounter<NemoFieldReader> {
 public:
  static const std::string classname() {return "orcamodel::NemoFieldReader";}

  explicit NemoFieldReader(eckit::PathName& filename);

  std::vector<atlas::PointXY> read_locs();
  size_t read_dim_size(const std::string& name);
  void read_datetimes();
  size_t get_nearest_datetime_index(const util::DateTime& datetime);
  template<class T> T read_fillvalue(const std::string& name);
  std::vector<double> read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx);
  void read_surf_var(const std::string& varname, const size_t t_indx,
      atlas::array::ArrayView<double, 2>& field_view);
  void read_volume_var(const std::string& varname, const size_t t_indx,
      atlas::array::ArrayView<double, 2>& field_view);
  void read_vertical_var(const std::string& varname,
      atlas::array::ArrayView<double, 2>& field_view);
  void read_vertical_var(const std::string& varname, const atlas::Mesh& mesh,
      atlas::array::ArrayView<double, 2>& field_view);

  void read_surf_var(const std::string& varname, const atlas::Mesh& mesh,
      const size_t t_indx, atlas::array::ArrayView<double, 2>& field_view);
  void read_volume_var(const std::string& varname,
     const atlas::Mesh& mesh, const size_t t_indx,
     atlas::array::ArrayView<double, 2>& field_view);

 private:
  NemoFieldReader() : ncFile() {}
  std::unique_ptr<netCDF::NcFile> ncFile;
  std::vector<util::DateTime> datetimes_;
  std::string time_dimvar_name_;
  std::string z_dimvar_name_;
};
}  // namespace orcamodel

//------------------------------------------------------------------------------
