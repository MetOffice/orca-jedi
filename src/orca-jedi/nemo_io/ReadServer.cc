/*
 * (C) British Crown Copyright 2026 Met Office
 */



#include "orca-jedi/nemo_io/ReadServer.h"

#include <string>
#include <vector>
#include <memory>

#include "oops/util/Logger.h"
#include "atlas/parallel/omp/omp.h"
#include "eckit/exception/Exceptions.h"

#include "orca-jedi/nemo_io/ParallelNetCDFBackend.h"

namespace orcamodel {

ReadServer::ReadServer(std::shared_ptr<eckit::Timer> eckit_timer,
  const eckit::PathName& file_path, const atlas::Mesh& mesh,
  bool parallel_io, size_t n_io_ranks) :
  mesh_(mesh),
  eckit_timer_(eckit_timer),
  parallel_(parallel_io) {
  buffer_indices_ = AtlasIndexToBufferIndexCreator::create_unique(
          mesh.grid().type(), mesh);

  // The root rank always keeps a serial reader for cheap metadata reads
  // (datetime index, fill value) and small 1-D vertical variables. Only the
  // large gridded field reads are routed through the parallel pool.
  if (myrank == mpiroot) {
    reader_ = std::make_unique<NemoFieldReader>(file_path);
  }

  if (parallel_) {
    this->setup_parallel(file_path, n_io_ranks);
  }
}

/// \brief Read in a 2D horizontal slice of variable data on the root processor only
/// \param var_name
/// \param t_index Index of the time slice.
/// \param z_index Index of the vertical slice.
/// \param buffer Vector to store the data.
template<class T> void ReadServer::read_var_on_root(const std::string& var_name,
              const size_t t_index,
              const size_t z_index,
              std::vector<T>& buffer) const {
  oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::read_var_on_root "
    << var_name << std::endl;
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();
  if (myrank == mpiroot) {
    buffer = reader_->read_var_slice<T>(var_name, t_index, z_index);
  } else {
    buffer.resize(size);
  }
}
template void ReadServer::read_var_on_root<double>(const std::string& var_name,
              const size_t t_index,
              const size_t z_index,
              std::vector<double>& buffer) const;
template void ReadServer::read_var_on_root<float>(const std::string& var_name,
              const size_t t_index,
              const size_t z_index,
              std::vector<float>& buffer) const;

/// \brief Read 1D vertical variable data on the root processor only
/// \param var_name NetCDF name of the vertical variable.
/// \param n_levels Number of levels to read from the file
/// \param buffer Vector to store the data.
template<class T> void ReadServer::read_vertical_var_on_root(const std::string& var_name,
              const size_t n_levels,
              std::vector<T>& buffer) const {
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();
  if (myrank == mpiroot) {
    buffer = reader_->read_vertical_var<T>(var_name, n_levels);
  } else {
    buffer.resize(size);
  }
}
template void ReadServer::read_vertical_var_on_root<double>(const std::string& var_name,
              const size_t n_levels,
              std::vector<double>& buffer) const;
template void ReadServer::read_vertical_var_on_root<float>(const std::string& var_name,
              const size_t n_levels,
              std::vector<float>& buffer) const;

/// \brief Move data from a buffer into an atlas arrayView.
/// \param buffer Vector of data to read
/// \param z_index Index of the vertical slice.
/// \param field_view View into the atlas field to store the data.
template<class T> void ReadServer::fill_field(const std::vector<T>& buffer,
              const size_t z_index,
      atlas::array::ArrayView<T, 2>& field_view) const {
    oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::fill_field" << std::endl;
    auto ghost = atlas::array::make_view<int32_t, 1>(this->mesh_.nodes().ghost());
    const size_t num_nodes = field_view.shape(0);
    // "ReadServer buffer size does not equal the number of horizontal nodes in the field_view"
    ASSERT(num_nodes <= buffer.size());
    atlas_omp_parallel_for(size_t inode = 0; inode < num_nodes; ++inode) {
      if (ghost(inode)) continue;
      const int64_t ibuf = (*buffer_indices_)(inode);
      field_view(inode, z_index) = buffer[ibuf];
    }
}

template void ReadServer::fill_field<double>(
      const std::vector<double>& buffer,
      const size_t z_index,
      atlas::array::ArrayView<double, 2>& field_view) const;
template void ReadServer::fill_field<float>(
      const std::vector<float>& buffer,
      const size_t z_index,
      atlas::array::ArrayView<float, 2>& field_view) const;

/// \brief Move vertical data from a buffer into an atlas arrayView.
/// \param buffer Vector of data to read
/// \param field_view View into the atlas field to store the data.
template<class T> void ReadServer::fill_vertical_field(const std::vector<T>& buffer,
      atlas::array::ArrayView<T, 2>& field_view) const {
    oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::fill_vertical_field "<< std::endl;
    auto ghost = atlas::array::make_view<int32_t, 1>(this->mesh_.nodes().ghost());
    const size_t num_nodes = field_view.shape(0);
    const size_t num_levels = field_view.shape(1);
    // "ReadServer buffer size does not equal the number of levels in the field_view"
    ASSERT(num_levels <= buffer.size());

    // even for 1D depths, store the data in an atlas 3D field - inefficient but flexible
    // NOTE: Use thread-private levels data to reduce OMP cache misses
    atlas_omp_parallel {
      const std::vector<T> buffer_TP = buffer;
      atlas_omp_for(size_t inode = 0; inode < num_nodes; ++inode) {
        for (size_t k = 0; k < num_levels; ++k) {
          if (ghost(inode)) continue;
          field_view(inode, k) = buffer_TP[k];
        }
      }
    }
}

template void ReadServer::fill_vertical_field<double>(
      const std::vector<double>& buffer,
      atlas::array::ArrayView<double, 2>& field_view) const;
template void ReadServer::fill_vertical_field<float>(
      const std::vector<float>& buffer,
      atlas::array::ArrayView<float, 2>& field_view) const;

/// \brief Read a NetCDF variable into an atlas field.
/// \param var_name The netCDF name of the variable to read.
/// \param t_index The time index for the data to read.
/// \param field_view View into the atlas field to store the data.
template<class T> void ReadServer::read_var(const std::string& var_name,
    const size_t t_index,
    atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::read_var "
    << var_name << std::endl;

  if (parallel_) {
    this->read_var_parallel<T>(var_name, t_index, field_view);
    return;
  }

  size_t n_levels = field_view.shape(1);
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();

  std::vector<T> buffer;
  // For each level
  for (size_t iLev = 0; iLev < n_levels; iLev++) {
    // read the data for that level onto the root processor
    this->read_var_on_root<T>(var_name, t_index, iLev, buffer);

    ASSERT(buffer.size() == size);
    // mpi distribute that data out to all processors
    atlas::mpi::comm().broadcast(&buffer.front(), size, mpiroot);

    // each processor fills out its field_view from the buffer
    this->fill_field<T>(buffer, iLev, field_view);
    log_status();
  }
}

template void ReadServer::read_var<double>(
    const std::string& var_name,
    const size_t t_index,
    atlas::array::ArrayView<double, 2>& field_view);
template void ReadServer::read_var<float>(
    const std::string& var_name,
    const size_t t_index,
    atlas::array::ArrayView<float, 2>& field_view);

// -----------------------------------------------------------------------------
// Parallel (I/O pool + redistributor + parallel-netCDF backend) read path.
// -----------------------------------------------------------------------------

void ReadServer::setup_parallel(const eckit::PathName& file_path,
    size_t n_io_ranks) {
  oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::setup_parallel" << std::endl;
  const size_t nx = buffer_indices_->nx();
  const size_t ny = buffer_indices_->ny();

  // Every local node (including ghost nodes) participates so the scatter fills
  // the ORCA halo cells that are only reachable through ghost nodes.
  std::vector<size_t> node_index;
  std::vector<size_t> global_index;
  node_index.reserve(mesh_.nodes().size());
  global_index.reserve(mesh_.nodes().size());
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    node_index.emplace_back(static_cast<size_t>(i_node));
    global_index.emplace_back((*buffer_indices_)(i_node));
  }

  const eckit::mpi::Comm& comm = atlas::mpi::comm();
  const size_t requested = n_io_ranks == 0 ? comm.size() : n_io_ranks;
  pool_ = std::make_unique<IoPool>(comm, requested, nx, ny, "orca_read_pool");
  redist_ = std::make_unique<Redistributor>(comm, *pool_, node_index, global_index);

  if (pool_->is_io_rank()) {
    backend_ = std::make_unique<ParallelNetCDFReadBackend>(pool_->io_comm(), file_path);
  }
}

/// \brief Read a variable one level at a time through the I/O pool and scatter
///        each level onto the field view. Unlike the root-broadcast path this
///        also fills ghost nodes (a superset of what the caller needs; the
///        subsequent halo exchange is a no-op on those points).
template<class T> void ReadServer::read_var_parallel(const std::string& var_name,
    const size_t t_index, atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::read_var_parallel "
    << var_name << std::endl;
  const bool is_io = pool_->is_io_rank();
  const Hyperslab slab = is_io ? pool_->decomposition().slab(pool_->io_rank())
                               : Hyperslab{};
  const size_t n_levels = field_view.shape(1);
  const size_t num_nodes = field_view.shape(0);

  std::vector<T> local(num_nodes);
  for (size_t iLev = 0; iLev < n_levels; ++iLev) {
    std::vector<T> slab_data;
    if (is_io) {
      backend_->read_slab(var_name, t_index, iLev, slab, slab_data);
    }
    redist_->from_io<T>(slab_data, local);
    for (size_t inode = 0; inode < num_nodes; ++inode) {
      field_view(inode, iLev) = local[inode];
    }
    log_status();
  }
}
template void ReadServer::read_var_parallel<double>(const std::string& var_name,
    const size_t t_index, atlas::array::ArrayView<double, 2>& field_view);
template void ReadServer::read_var_parallel<float>(const std::string& var_name,
    const size_t t_index, atlas::array::ArrayView<float, 2>& field_view);

/// \brief Read a vertical variable into an atlas field.
/// \param var_name The netCDF name of the variable to read.
/// \param field_view View into the atlas field to store the data.
template<class T> void ReadServer::read_vertical_var(const std::string& var_name,
    atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::ReadServer::read_vertical_var "
    << var_name << std::endl;

  size_t n_levels = field_view.shape(1);

  std::vector<T> buffer;

  // read the data onto the root processor
  this->read_vertical_var_on_root<T>(var_name, n_levels, buffer);

  // mpi distribute that data out to all processors
  atlas::mpi::comm().broadcast(&buffer.front(), n_levels, mpiroot);

  // each processor fills out its field_view from the buffer
  this->fill_vertical_field<T>(buffer, field_view);
  log_status();
}

template void ReadServer::read_vertical_var<double>(
    const std::string& var_name,
    atlas::array::ArrayView<double, 2>& field_view);
template void ReadServer::read_vertical_var<float>(
    const std::string& var_name,
    atlas::array::ArrayView<float, 2>& field_view);

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
