/*
 * (C) British Crown Copyright 2024 Met Office
 */



#include "orca-jedi/nemo_io/ReadServer.h"
#include <string>
#include <vector>
#include <memory>

#include "eckit/exception/Exceptions.h"

namespace orcamodel {

ReadServer::ReadServer(const eckit::PathName& file_path, const atlas::Mesh& mesh) :
  mesh_(mesh),
  index_glbarray_(atlas::OrcaGrid(mesh.grid())) {
  if (myrank == mpiroot) {
    reader_ = std::make_unique<NemoFieldReader>(file_path);
  }
}

/// \brief Read in a 2D horizontal slice of variable data on the root processor only
/// \param var_name
/// \param t_index Index of the time slice.
/// \param z_index Index of the vertical slice.
/// \param buffer Vector to store the data.
void ReadServer::read_var_on_root(const std::string& var_name,
              const size_t t_index,
              const size_t z_index,
              std::vector<double>& buffer) const {
  size_t size = index_glbarray_.nx_halo_WE * index_glbarray_.ny_halo_NS;
  if (myrank == mpiroot) {
    buffer = reader_->read_var_slice(var_name, t_index, z_index);
  } else {
    buffer.resize(size);
  }
}

/// \brief Read 1D vertical variable data on the root processor only
/// \param var_name NetCDF name of the vertical variable.
/// \param n_levels Number of levels to read from the file
/// \param buffer Vector to store the data.
void ReadServer::read_vertical_var_on_root(const std::string& var_name,
              const size_t n_levels,
              std::vector<double>& buffer) const {
  size_t size = index_glbarray_.nx_halo_WE * index_glbarray_.ny_halo_NS;
  if (myrank == mpiroot) {
    buffer = reader_->read_vertical_var(var_name, n_levels);
  } else {
    buffer.resize(size);
  }
}

/// \brief Move data from a buffer into an atlas arrayView.
/// \param buffer Vector of data to read
/// \param z_index Index of the vertical slice.
/// \param field_view View into the atlas field to store the data.
void ReadServer::fill_field(const std::vector<double>& buffer,
              const size_t z_index,
      atlas::array::ArrayView<double, 2>& field_view) const {
    auto ghost = atlas::array::make_view<int32_t, 1>(this->mesh_.nodes().ghost());
    auto ij = atlas::array::make_view<int32_t, 2>(this->mesh_.nodes().field("ij"));
    const size_t numNodes = field_view.shape(0);
    // "ReadServer buffer size does not equal the number of horizontal nodes in the field_view"
    assert(numNodes <= buffer.size());
    for (size_t inode = 0; inode < numNodes; ++inode) {
      if (ghost(inode)) continue;
      const int64_t ibuf = index_glbarray_(ij(inode, 0), ij(inode, 1));
      field_view(inode, z_index) = buffer[ibuf];
    }
}

/// \brief Move vertical data from a buffer into an atlas arrayView.
/// \param buffer Vector of data to read
/// \param field_view View into the atlas field to store the data.
void ReadServer::fill_vertical_field(const std::vector<double>& buffer,
      atlas::array::ArrayView<double, 2>& field_view) const {
    auto ghost = atlas::array::make_view<int32_t, 1>(this->mesh_.nodes().ghost());
    const size_t num_nodes = field_view.shape(0);
    const size_t num_levels = field_view.shape(1);
    // "ReadServer buffer size does not equal the number of levels in the field_view"
    assert(num_levels <= buffer.size());

    // even for 1D depths, store the data in an atlas 3D field - inefficient but flexible
    for (size_t inode = 0; inode < num_nodes; ++inode) {
      for (size_t k = 0; k < num_levels; ++k) {
        if (ghost(inode)) continue;
        field_view(inode, k) = buffer[k];
      }
    }
}

/// \brief Read a NetCDF variable into an atlas field.
/// \param var_name The netCDF name of the variable to read.
/// \param t_index The time index for the data to read.
/// \param field_view View into the atlas field to store the data.
void ReadServer::read_var(const std::string& var_name,
    const size_t t_index,
    atlas::array::ArrayView<double, 2>& field_view) {

  size_t n_levels = field_view.shape(1);
  size_t size = index_glbarray_.nx_halo_WE * index_glbarray_.ny_halo_NS;

  std::vector<double> buffer;
  // For each level
  for (size_t iLev = 0; iLev < n_levels; iLev++) {
    // read the data for that level onto the root processor
    this->read_var_on_root(var_name, t_index, iLev, buffer);

    // mpi distribute that data out to all processors
    atlas::mpi::comm().broadcast(&buffer.front(), size, mpiroot);

    // each processor fills out its field_view from the buffer
    this->fill_field(buffer, iLev, field_view);
  }
}

/// \brief Read a vertical variable into an atlas field.
/// \param var_name The netCDF name of the variable to read.
/// \param field_view View into the atlas field to store the data.
void ReadServer::read_vertical_var(const std::string& var_name,
    atlas::array::ArrayView<double, 2>& field_view) {

  size_t n_levels = field_view.shape(1);
  size_t size = index_glbarray_.nx_halo_WE * index_glbarray_.ny_halo_NS;

  std::vector<double> buffer;

  // read the data onto the root processor
  this->read_vertical_var_on_root(var_name, n_levels, buffer);

  // mpi distribute that data out to all processors
  atlas::mpi::comm().broadcast(&buffer.front(), size, mpiroot);

  // each processor fills out its field_view from the buffer
  this->fill_vertical_field(buffer, field_view);
}

/// \brief Find the nearest datetime index to a datetime on the MPI root only.
/// \param datetime Search for the index of the time slice in the file nearest this datetime.
/// \return The index of the nearest time slice in the file.
size_t ReadServer::get_nearest_datetime_index(const util::DateTime& datetime) const {
  size_t t_index;

  if (myrank == mpiroot) {
    t_index = reader_->get_nearest_datetime_index(datetime);
  }

  // mpi distribute that data out to all processors
  atlas::mpi::comm().broadcast(t_index, mpiroot);

  return t_index;
}

/// \brief Read the _FillValue for a variable, defaulting to the minimum value
/// for the datatype. Read on the MPI root process only.
/// \param name Name of the netCDF variable containing the _FillValue attribute to retrieve.
/// \return The fill value for this netCDF variable.
template<class T> T ReadServer::read_fillvalue(const std::string& name) const {
  T fillvalue;

  if (myrank == mpiroot) {
    fillvalue = reader_->read_fillvalue<T>(name);
  }

  // mpi distribute that data out to all processors
  atlas::mpi::comm().broadcast(fillvalue, mpiroot);

  return fillvalue;
}
template int ReadServer::read_fillvalue<int>(const std::string& name) const;
template float ReadServer::read_fillvalue<float>(const std::string& name) const;
template double ReadServer::read_fillvalue<double>(const std::string& name) const;
}  // namespace orcamodel
