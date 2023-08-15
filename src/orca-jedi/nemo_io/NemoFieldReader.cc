/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include <netcdf>
// Using Lynton Appel's netcdf-cxx4 from
// https://github.com/Unidata/netcdf-cxx4

#include <algorithm>
#include <sstream>
#include <limits>
#include <map>

#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

#include "atlas/field.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

namespace {
struct IndexGlbArray {
    int32_t ix_glb_max;
    int32_t iy_glb_max;
    int32_t glbarray_offset;
    int32_t glbarray_jstride;
    int32_t nx_halo_WE;
    int32_t ny_halo_NS;

    explicit IndexGlbArray(const atlas::OrcaGrid& orcaGrid) {
        iy_glb_max = orcaGrid.ny() + orcaGrid.haloNorth() - 1;
        ix_glb_max = orcaGrid.nx() + orcaGrid.haloEast() - 1;

        nx_halo_WE = orcaGrid.nx() + orcaGrid.haloEast() + orcaGrid.haloWest();
        ny_halo_NS = orcaGrid.ny() + orcaGrid.haloNorth()
          + orcaGrid.haloSouth();

        // vector of local indices: necessary for remote indices of ghost nodes
        int iy_glb_min = -orcaGrid.haloSouth();
        int ix_glb_min = -orcaGrid.haloWest();
        glbarray_offset  = -(nx_halo_WE * iy_glb_min) - ix_glb_min;
        glbarray_jstride = nx_halo_WE;
    }

    int64_t operator()(int i, int j) {
        ATLAS_ASSERT(i <= ix_glb_max,
            std::to_string(i) + " > " + std::to_string(ix_glb_max));
        ATLAS_ASSERT(j <= iy_glb_max,
            std::to_string(j) + " > " + std::to_string(iy_glb_max));
        return glbarray_offset + j * glbarray_jstride + i;
    }
};

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
}  // namespace

NemoFieldReader::NemoFieldReader(eckit::PathName& filename)
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

size_t NemoFieldReader::read_dim_size(const std::string& name) {
  auto dim = ncFile->getDim(name);
  if (dim.isNull()) {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_dim_size Dimension '"
               << name << "' is not present in NetCDF file" << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  oops::Log::debug() << "orcamodel::NemoFieldReader:: group name "
                     << ncFile->getName(true) << " dim name: "  << name
                     << std::endl;

  return dim.getSize();
}

/// \brief retrieve the datetime for each time index in file.
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
template<class T> T NemoFieldReader::read_fillvalue(const std::string& name) {
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
template int NemoFieldReader::read_fillvalue<int>(const std::string& name);
template float NemoFieldReader::read_fillvalue<float>(const std::string& name);
template double NemoFieldReader::read_fillvalue<double>(
    const std::string& name);

/// \brief get the time dimension index corresponding to the nearest datetime
/// to a target datetime.
size_t NemoFieldReader::get_nearest_datetime_index(
    const util::DateTime& tgt_datetime) {
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

std::vector<atlas::PointXY> NemoFieldReader::read_locs() {
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

    nc_var_lat.getVar({0, 0}, {ny, nx}, lats.data());

    varname = "nav_lon";
    netCDF::NcVar nc_var_lon = ncFile->getVar(varname);
    if (nc_var_lon.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_locs ncVar '" << varname
                 <<"' is not present in NetCDF file.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> lons(nx*ny);

    nc_var_lon.getVar({0, 0}, {ny, nx}, lons.data());

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

std::vector<double> NemoFieldReader::read_var_slice(const std::string& varname,
      const size_t t_indx, const size_t z_indx) {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_var_slice"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> var_data(nx*ny);

    size_t n_dims = nc_var.getDimCount();
    if (n_dims == 4) {
      nc_var.getVar({t_indx, z_indx, 0, 0}, {1, 1, ny, nx}, var_data.data());
    } else if (n_dims == 3) {
      nc_var.getVar({t_indx, 0, 0}, {1, ny, nx}, var_data.data());
    } else if (n_dims == 2) {
      nc_var.getVar({0, 0}, {ny, nx}, var_data.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    return var_data;
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_surf_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

  void NemoFieldReader::read_volume_var(const std::string& varname,
    const size_t t_indx, atlas::array::ArrayView<double, 2>& field_view) {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_volume_var"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");
    size_t nz = read_dim_size(z_dimvar_name_);
    size_t nlevels = field_view.shape(1);

    if (field_view.shape(0) != nx*ny) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var field_view 1st"
                 << " dimension does not match horizontal dimensions"
                 << " for varname " << varname;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    if (nlevels > nz) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var field_view 2nd"
                 << " dimension " << nlevels << " is larger than NetCDF file"
                 << " z dimension " << nz << " for varname " << varname;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> buffer(nx*ny*nlevels);

    size_t n_dims = nc_var.getDimCount();
    std::string first_dim_name = nc_var.getDim(0).getName();
    if (n_dims == 4) {
      nc_var.getVar({t_indx, 0, 0, 0}, {1, nlevels, ny, nx}, buffer.data());
    } else if (n_dims == 3 && first_dim_name == z_dimvar_name_) {
      nc_var.getVar({0, 0, 0}, {nlevels, ny, nx}, buffer.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    // in atlas fields the levels indices change the fastest, so we need to
    // swap the indexing order from the netCDF data.
    for (int n = 0; n < nx*ny; ++n) {
      for (int k = 0; k < nlevels; ++k) {
        field_view(n, k) = buffer[k*nx*ny + n];
      }
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_volume_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

void NemoFieldReader::read_vertical_var(const std::string& varname,
    atlas::array::ArrayView<double, 2>& field_view) {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_vertical_var"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");
    size_t nz = read_dim_size(z_dimvar_name_);
    size_t nlevels = field_view.shape(1);

    if (nlevels > nz) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var field_view "
                 << "2nd dimension is larger than NetCDF file z dimension ";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> buffer(nlevels);

    size_t n_dims = nc_var.getDimCount();
    std::string first_dim_name = nc_var.getDim(0).getName();
    if (n_dims == 1 && first_dim_name == z_dimvar_name_) {
      nc_var.getVar({0}, {nlevels}, buffer.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    // Store the data in an atlas 3D field - inefficient but flexible
    for (size_t inode = 0; inode < field_view.shape(0); ++inode) {
      for (int k = 0; k < nlevels; ++k) {
        field_view(inode, k) = buffer[k];
      }
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_vertical_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

void NemoFieldReader::read_vertical_var(const std::string& varname,
    const atlas::Mesh& mesh, atlas::array::ArrayView<double, 2>& field_view) {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_vertical_var"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");
    size_t nz = read_dim_size(z_dimvar_name_);
    size_t nlevels = field_view.shape(1);

    if (nlevels > nz) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var field_view "
                 << "2nd dimension is larger than NetCDF file z dimension ";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> buffer(nlevels);

    size_t n_dims = nc_var.getDimCount();
    std::string first_dim_name = nc_var.getDim(0).getName();
    if (n_dims == 1 && first_dim_name == z_dimvar_name_) {
      nc_var.getVar({0}, {nlevels}, buffer.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_vertical_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    // Store the data in an atlas 3D field - inefficient but flexible
    for (size_t inode = 0; inode < field_view.shape(0); ++inode) {
      for (int k = 0; k < nlevels; ++k) {
        if (ghost(inode)) continue;
        field_view(inode, k) = buffer[k];
      }
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_vertical_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

void NemoFieldReader::read_surf_var(const std::string& varname,
    const size_t t_indx, atlas::array::ArrayView<double, 2>& field_view) {
  oops::Log::trace() << "orcamodel::NemoFieldReader::read_surf_var"
                     << std::endl;
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");

    if (field_view.size() != nx*ny) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var field_view "
                 << "dimensions  do not match dimensions in netCDF file ";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    size_t n_dims = nc_var.getDimCount();
    if (n_dims == 4) {
      nc_var.getVar({t_indx, 0, 0, 0}, {1, 1, ny, nx}, field_view.data());
    } else if (n_dims == 3) {
      nc_var.getVar({t_indx, 0, 0}, {1, ny, nx}, field_view.data());
    } else if (n_dims == 2) {
      nc_var.getVar({0, 0}, {ny, nx}, field_view.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_surf_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}

void NemoFieldReader::read_surf_var(const std::string& varname,
    const atlas::Mesh& mesh, const size_t t_indx,
    atlas::array::ArrayView<double, 2>& field_view) {
    oops::Log::trace() << "orcamodel::NemoFieldReader::read_surf_var"
                       << std::endl;
    try {
        size_t nx = read_dim_size("x");
        size_t ny = read_dim_size("y");

        const atlas::OrcaGrid orcaGrid = atlas::OrcaGrid(mesh.grid());
        if (!orcaGrid) {
          std::ostringstream err_stream;
          err_stream << "orcamodel::NemoFieldReader::read_surf_var"
                     << " only reads ORCA grid data. " << mesh.grid().name()
                     << " is not an ORCA grid." << std::endl;
          throw eckit::BadValue(err_stream.str(), Here());
        }

        IndexGlbArray index_glbarray(orcaGrid);
        if (index_glbarray.nx_halo_WE != nx
            || index_glbarray.ny_halo_NS != ny) {
            std::ostringstream err_stream;
            err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar"
                << " grid dimensions don't match file dimensions grid: ("
                << index_glbarray.nx_halo_WE << ", "
                << index_glbarray.ny_halo_NS << ") file: " << nx << ", "
                << ny << ")" << std::endl;
            throw eckit::BadValue(err_stream.str());
        }

        auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
        auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

        if (mesh.nodes().size() > nx * ny) {
            std::ostringstream err_stream;
            err_stream << "orcamodel::NemoFieldReader::read_surf_var ncVar"
                << " number of mesh nodes " << "larger than the netCDF file"
                << " dimensions, " << nx * ny
                << " double check the grids match." << std::endl;
            throw eckit::BadValue(err_stream.str());
        }

        netCDF::NcVar nc_var = ncFile->getVar(varname);
        if (nc_var.isNull()) {
            throw eckit::UserError(
                std::string("orcamodel::NemoFieldReader::read_surf_var ncVar ")
                + varname + " is not present in NetCDF file", Here());
        }

        std::vector<double> var_data(nx * ny);
        size_t n_dims = nc_var.getDimCount();
        if (n_dims == 4) {
            nc_var.getVar({t_indx, 0, 0, 0}, {1, 1, ny, nx}, var_data.data());
        } else if (n_dims == 3) {
           nc_var.getVar({t_indx, 0, 0}, {1, ny, nx}, var_data.data());
        } else if (n_dims == 2) {
           nc_var.getVar({0, 0}, {ny, nx}, var_data.data());
        } else {
            throw eckit::UserError(
                std::string("orcamodel::NemoFieldReader::read_surf_var ")
                + "ncVar " + varname + " has an unreadable number of "
                + "dimensions: " + std::to_string(n_dims), Here());
        }

        for (size_t inode = 0; inode < field_view.size(); ++inode) {
          if (ghost(inode)) continue;
          double data = var_data[index_glbarray(ij(inode, 0), ij(inode, 1))];
          field_view(inode, 0) = data;
        }
    }
    catch (netCDF::exceptions::NcException& e) {
        throw eckit::FailedLibraryCall("NetCDF",
            "orcamodel::NemoFieldReader::read_surf_var", e.what(), Here());
    }
}

void NemoFieldReader::read_volume_var(const std::string& varname,
    const atlas::Mesh& mesh, const size_t t_indx,
    atlas::array::ArrayView<double, 2>& field_view) {
  try {
    size_t nx = read_dim_size("x");
    size_t ny = read_dim_size("y");
    size_t nz = read_dim_size(z_dimvar_name_);
    size_t nlevels = field_view.shape(1);

    const atlas::OrcaGrid orcaGrid = atlas::OrcaGrid(mesh.grid());
    if (!orcaGrid) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var"
                 << " only reads ORCA grid data. " << mesh.grid().name()
                 << " is not an ORCA grid." << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    IndexGlbArray index_glbarray(orcaGrid);
    if (index_glbarray.nx_halo_WE != nx || index_glbarray.ny_halo_NS != ny) {
        std::ostringstream err_stream;
        err_stream << "orcamodel::NemoFieldReader::read_volume_var ncVar"
            << " grid dimensions don't match file dimensions grid: ("
            << index_glbarray.nx_halo_WE << ", "
            << index_glbarray.ny_halo_NS << ") file: " << nx << ", "
            << ny << ")" << std::endl;
        throw eckit::BadValue(err_stream.str());
    }

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    if (mesh.nodes().size() > nx * ny) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var field_view 1st"
                 << " dimension does not match horizontal dimensions"
                 << " for varname " << varname;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    if (nlevels > nz) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var field_view 2nd"
                 << " dimension " << nlevels << " is larger than NetCDF file"
                 << " z dimension " << nz << " for varname " << varname;
      throw eckit::BadValue(err_stream.str(), Here());
    }

    netCDF::NcVar nc_var = ncFile->getVar(varname);
    if (nc_var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var ncVar '"
                 << varname << "' is not present in NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    std::vector<double> buffer(nx*ny*nlevels);

    size_t n_dims = nc_var.getDimCount();
    std::string first_dim_name = nc_var.getDim(0).getName();
    if (n_dims == 4) {
      nc_var.getVar({t_indx, 0, 0, 0}, {1, nlevels, ny, nx}, buffer.data());
    } else if (n_dims == 3 && first_dim_name == z_dimvar_name_) {
      nc_var.getVar({0, 0, 0}, {nlevels, ny, nx}, buffer.data());
    } else {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFieldReader::read_volume_var ncVar '"
                 << varname << "' has " << n_dims << " dimensions.";
      throw eckit::BadValue(err_stream.str(), Here());
    }

    // in atlas fields the levels indices change the fastest, so we need to
    // swap the indexing order from the netCDF data.
    const size_t numNodes = field_view.shape(0);
    atlas_omp_for(int k = 0; k < nlevels; ++k) {
      for (size_t inode = 0; inode < numNodes; ++inode) {
        if (ghost(inode)) continue;
        field_view(inode, k) =
          buffer[k*nx*ny + index_glbarray(ij(inode, 0), ij(inode, 1))];
      }
    }
  } catch(netCDF::exceptions::NcException& e)
  {
    std::ostringstream err_stream;
    err_stream << "orcamodel::NemoFieldReader::read_volume_var varname: "
               << varname << " NetCDF exception: " << std::endl << e.what();
    throw eckit::ReadError(err_stream.str(), Here());
  }
}
}  // namespace orcamodel
