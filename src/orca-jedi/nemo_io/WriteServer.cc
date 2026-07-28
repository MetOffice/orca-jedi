/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/WriteServer.h"

#include<algorithm>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "atlas-orca/grid/OrcaGrid.h"

#include "orca-jedi/nemo_io/ParallelNetCDFBackend.h"


namespace orcamodel {

/// \brief  Write a 3D field at a given time index to a file on the root process.
/// \param var_name   Name of the variable in the file.
/// \param t_index    Time index in the file to write.
/// \param buffer     vector of data to write.
template<class T> void WriteServer::write_vol_var_on_root(const std::string& var_name,
    const size_t t_index,
    const std::vector<T>& buffer) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_vol_var_on_root "
                     << std::endl;
  if (myrank == mpiroot) {
    writer_->write_vol_var<T>(var_name, buffer, t_index);
  }
}
template void WriteServer::write_vol_var_on_root<double>(const std::string& var_name,
    const size_t t_index, const std::vector<double>& buffer);
template void WriteServer::write_vol_var_on_root<float>(const std::string& var_name,
    const size_t t_index, const std::vector<float>& buffer);

/// \brief  Write a 2D field at a given time index to a file on the root process.
/// \param var_name   Name of the variable in the file.
/// \param t_index    Time index in the file to write.
/// \param buffer     vector of data to write.
template<class T> void WriteServer::write_surf_var_on_root(const std::string& var_name,
    const size_t t_index,
    const std::vector<T>& buffer) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_surf_var_on_root "
                     << std::endl;
  if (myrank == mpiroot) {
    writer_->write_surf_var<T>(var_name, buffer, t_index);
  }
}
template void WriteServer::write_surf_var_on_root<double>(const std::string& var_name,
    const size_t t_index, const std::vector<double>& buffer);
template void WriteServer::write_surf_var_on_root<float>(const std::string& var_name,
    const size_t t_index, const std::vector<float>& buffer);

/// \brief sort the data into the correct order for the NetCDF file
/// \param  buffer        .
/// \return sorted_buffer vector of data to write.
template<class T> std::vector<T> WriteServer::sort_buffer(
  const std::vector<T> & buffer) const {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::sort_buffer "
                     << std::endl;
  std::vector<T> sorted_buffer(buffer_indices_->nx() * buffer_indices_->ny());
  for (size_t i_buf = 0; i_buf < buffer.size(); ++i_buf) {
    const size_t source_index = unsorted_buffer_indices_[i_buf];
    ASSERT(source_index < sorted_buffer.size());
    sorted_buffer[source_index] = buffer[i_buf];
  }
  return sorted_buffer;
}
template std::vector<double> WriteServer::sort_buffer<double>(
    const std::vector<double> & buffer) const;
template std::vector<float> WriteServer::sort_buffer<float>(
    const std::vector<float> & buffer) const;

/// \brief Gather field_view data from all processes to a buffer on the server root process
/// \param field_view   Field data to write to the file.
/// \param i_level      level index for the data.
/// \param missingValue missing value object for matching masked data.
/// \return buffer      vector of data to write.
template<class T> std::vector<T> WriteServer::gather_field_on_root(
  const atlas::array::ArrayView<T, 2>& field_view, const size_t i_level,
  const atlas::field::MissingValue& missingValue) const {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::gather_field_on_root "
                     << std::endl;

  std::vector<T> local_buffer;
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();
  const bool has_mv = static_cast<bool>(missingValue);

  // handle the serial data case
  if (is_serial_) {
    if (myrank == mpiroot) {
      for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
        if (has_mv && missingValue(field_view(i_node, i_level))) {
          local_buffer.emplace_back(NemoFieldWriter::fillValue);
        } else {
          local_buffer.emplace_back(field_view(i_node, i_level));
        }
      }
      ASSERT(local_buffer.size() == size);
    }
    return local_buffer;
  }

  // handle the distributed data case
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    if (has_mv && missingValue(field_view(i_node, i_level))) {
      local_buffer.emplace_back(NemoFieldWriter::fillValue);
    } else {
      local_buffer.emplace_back(field_view(i_node, i_level));
    }
  }

  std::vector<T> unsorted_buffer(recvcnt_);
  atlas::mpi::comm().gatherv(&local_buffer.front(), local_buffer.size(),
                             &unsorted_buffer.front(), &recvcounts_.front(),
                             &recvdispls_.front(), mpiroot);

  if (myrank == mpiroot) {
    return this->sort_buffer<T>(unsorted_buffer);
  }

  return std::vector<T>();
}
template std::vector<double> WriteServer::gather_field_on_root<double>(
  const atlas::array::ArrayView<double, 2>& field_view, const size_t i_level,
  const atlas::field::MissingValue& missingValue) const;
template std::vector<float> WriteServer::gather_field_on_root<float>(
  const atlas::array::ArrayView<float, 2>& field_view, const size_t i_level,
  const atlas::field::MissingValue& missingValue) const;

/// \brief  Write a 3D field at a given time index to a file.
/// \param var_name     Name of the variable in the file.
/// \param t_index      Time index in the file to write.
/// \param missingValue missing value object for matching masked data.
/// \param field_view   Field data to write to the file.
template<class T> void WriteServer::write_vol_var(const std::string& var_name,
    const size_t t_index,
    const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_vol_var "
                     << var_name << std::endl;
  if (parallel_) {
    this->write_vol_var_parallel<T>(var_name, t_index, missingValue, field_view);
    return;
  }
  std::vector<T> buffer;
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();
  if (myrank == mpiroot) {
    buffer.resize(n_levels_*size);
  }
  // For each level
  for (size_t iLev = 0; iLev < n_levels_; iLev++) {
    // gather data for this level onto the root processor.
    const size_t start = iLev*size;
    const std::vector<T> local_buffer = this->gather_field_on_root<T>(field_view, iLev,
                                                                      missingValue);
    if (myrank == mpiroot) {
      ASSERT(local_buffer.size() == size);
      ASSERT(buffer.begin() + start + size <= buffer.end());
      std::copy(local_buffer.begin(), local_buffer.end(), buffer.begin()+start);
    }
    log_status();
  }
  if (myrank == mpiroot) {
    ASSERT(buffer.size() == n_levels_*size);
  }
  // write field from the buffer to variable
  this->write_vol_var_on_root<T>(var_name, t_index, buffer);
  log_status();
}
template void WriteServer::write_vol_var<double>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<double, 2>& field_view);
template void WriteServer::write_vol_var<float>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<float, 2>& field_view);

/// \brief  Write a 2D field at a given time index to a file.
/// \param var_name     Name of the variable in the file.
/// \param t_index      Time index in the file to write.
/// \param missingValue missing value object for matching masked data.
/// \param field_view   Field data to write to the file.
template<class T> void WriteServer::write_surf_var(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_surf_var "
                     << var_name << std::endl;
  if (parallel_) {
    this->write_surf_var_parallel<T>(var_name, t_index, missingValue, field_view);
    return;
  }
  // gather data for this level onto the root processor.
  const std::vector<T> buffer = this->gather_field_on_root<T>(field_view, 0, missingValue);
  log_status();
  // write field from the buffer to variable.
  this->write_surf_var_on_root<T>(var_name, t_index, buffer);
  log_status();
}
template void WriteServer::write_surf_var<double>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<double, 2>& field_view);
template void WriteServer::write_surf_var<float>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<float, 2>& field_view);

/// \brief  Write latitude and longitude dimensions.
void WriteServer::write_dimensions() {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_dimensions" << std::endl;
  if (parallel_) {
    this->write_dimensions_parallel();
    return;
  }
  atlas::array::ArrayView<double, 2> lonlat{atlas::array::make_view<double, 2>(
                                              mesh_.nodes().lonlat())};
  atlas::field::MissingValue missingValue(mesh_.nodes().lonlat());
  size_t size = buffer_indices_->nx() * buffer_indices_->ny();

  const std::vector<double> lon_buffer = this->gather_field_on_root(lonlat, 0, missingValue);
  const std::vector<double> lat_buffer = this->gather_field_on_root(lonlat, 1, missingValue);

  if (myrank == mpiroot) {
    ASSERT(lat_buffer.size() == size);
    ASSERT(lon_buffer.size() == size);
    writer_->write_dimensions(lat_buffer, lon_buffer);
  }
}

// -----------------------------------------------------------------------------
// Parallel (I/O pool + redistributor + parallel-netCDF backend) path.
// -----------------------------------------------------------------------------

/// \brief Build the I/O pool, the cached redistribution plan and, on I/O ranks,
///        the parallel-netCDF backend (which creates the file collectively).
void WriteServer::setup_parallel(const eckit::PathName& file_path,
    const std::vector<util::DateTime>& datetimes,
    const std::vector<double>& depths, size_t n_io_ranks) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::setup_parallel" << std::endl;
  const size_t nx = buffer_indices_->nx();
  const size_t ny = buffer_indices_->ny();

  // Every local node that maps onto the buffer (including in-buffer ghost
  // nodes) contributes its global buffer index so that the ORCA halo cells
  // reachable only through ghosts are covered. Out-of-buffer nodes (the
  // north-fold pivot ghost row) have no file cell and must not write one.
  std::vector<size_t> node_index;
  std::vector<size_t> global_index;
  node_index.reserve(mesh_.nodes().size());
  global_index.reserve(mesh_.nodes().size());
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    const size_t inode = static_cast<size_t>(i_node);
    if (!buffer_indices_->maps_to_buffer(inode)) continue;
    node_index.emplace_back(inode);
    global_index.emplace_back((*buffer_indices_)(inode));
  }

  const eckit::mpi::Comm& comm = atlas::mpi::comm();
  const size_t requested = n_io_ranks == 0 ? comm.size() : n_io_ranks;
  pool_ = std::make_unique<IoPool>(comm, requested, nx, ny, "orca_write_pool");
  redist_ = std::make_unique<Redistributor>(comm, *pool_, node_index, global_index);

  if (pool_->is_io_rank()) {
    backend_ = std::make_unique<ParallelNetCDFWriteBackend>(pool_->io_comm(),
        file_path, nx, ny, datetimes, depths);
  }
}

/// \brief Redistribute one horizontal level of a field onto this rank's I/O slab,
///        substituting the fill value for masked points (mirrors
///        gather_field_on_root). Returns an empty vector on non-I/O ranks.
template<class T> std::vector<T> WriteServer::field_to_slab(
    const atlas::array::ArrayView<T, 2>& field_view, const size_t i_level,
    const atlas::field::MissingValue& missingValue) const {
  const bool has_mv = static_cast<bool>(missingValue);
  std::vector<T> local(mesh_.nodes().size());
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    if (has_mv && missingValue(field_view(i_node, i_level))) {
      local[i_node] = NemoFieldWriter::fillValue;
    } else {
      local[i_node] = field_view(i_node, i_level);
    }
  }
  std::vector<T> slab;
  redist_->to_io<T>(local, slab);
  return slab;
}

/// \brief  Write latitude and longitude dimensions on the I/O pool.
void WriteServer::write_dimensions_parallel() {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_dimensions_parallel"
                     << std::endl;
  atlas::array::ArrayView<double, 2> lonlat{atlas::array::make_view<double, 2>(
                                              mesh_.nodes().lonlat())};
  atlas::field::MissingValue missingValue(mesh_.nodes().lonlat());

  // NEMO layout: nav_lon is component 0, nav_lat is component 1.
  const std::vector<double> lon_band = this->field_to_slab<double>(lonlat, 0, missingValue);
  const std::vector<double> lat_band = this->field_to_slab<double>(lonlat, 1, missingValue);

  if (pool_->is_io_rank()) {
    const Hyperslab slab = pool_->decomposition().slab(pool_->io_rank());
    backend_->write_dimensions(slab, lat_band, lon_band);
  }
}

/// \brief  Write a 3D field at a given time index through the I/O pool.
template<class T> void WriteServer::write_vol_var_parallel(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_vol_var_parallel "
                     << var_name << std::endl;
  const bool is_io = pool_->is_io_rank();
  const Hyperslab slab = is_io ? pool_->decomposition().slab(pool_->io_rank())
                               : Hyperslab{};
  // One level at a time keeps the slab buffer small; every I/O rank loops
  // identically so the collective backend writes stay in step.
  for (size_t iLev = 0; iLev < n_levels_; ++iLev) {
    const std::vector<T> slab_data = this->field_to_slab<T>(field_view, iLev, missingValue);
    if (is_io) {
      backend_->write_vol_slab(var_name, t_index, iLev, 1, slab, slab_data);
    }
    log_status();
  }
}
template void WriteServer::write_vol_var_parallel<double>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<double, 2>& field_view);
template void WriteServer::write_vol_var_parallel<float>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<float, 2>& field_view);

/// \brief  Write a 2D field at a given time index through the I/O pool.
template<class T> void WriteServer::write_surf_var_parallel(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<T, 2>& field_view) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::write_surf_var_parallel "
                     << var_name << std::endl;
  const std::vector<T> slab_data = this->field_to_slab<T>(field_view, 0, missingValue);
  if (pool_->is_io_rank()) {
    const Hyperslab slab = pool_->decomposition().slab(pool_->io_rank());
    backend_->write_surf_slab(var_name, t_index, slab, slab_data);
  }
  log_status();
}
template void WriteServer::write_surf_var_parallel<double>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<double, 2>& field_view);
template void WriteServer::write_surf_var_parallel<float>(const std::string& var_name,
    const size_t t_index, const atlas::field::MissingValue& missingValue,
    const atlas::array::ArrayView<float, 2>& field_view);

WriteServer::WriteServer(std::shared_ptr<eckit::Timer> eckit_timer,
    const eckit::PathName& file_path,
    const atlas::Mesh& mesh,
    const std::vector<util::DateTime> datetimes,
    const std::vector<double> depths,
    bool is_serial,
    bool parallel_io,
    size_t n_io_ranks) : mesh_(mesh), eckit_timer_(eckit_timer), is_serial_(is_serial),
  n_levels_(depths.size()), parallel_(parallel_io) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::WriteServer" << std::endl;
  buffer_indices_ = AtlasIndexToBufferIndexCreator::create_unique(
          mesh.grid().type(), mesh);

  if (parallel_) {
    this->setup_parallel(file_path, datetimes, depths, n_io_ranks);
    this->write_dimensions();
    return;
  }

  std::vector<size_t> local_buf_indices;

  if (is_serial_) {
    if (myrank == mpiroot) {
      writer_ = std::make_unique<NemoFieldWriter>(file_path, datetimes,
                                                  buffer_indices_->nx(), buffer_indices_->ny(),
                                                  depths);
    }
    this->write_dimensions();
    recvcnt_ = buffer_indices_->nx() * buffer_indices_->ny();
    return;
  }

  auto orcaGrid = atlas::OrcaGrid(mesh_.grid());
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    local_buf_indices.emplace_back((*buffer_indices_)(i_node));
  }

  // gather all remote buffer indices
  int sendcnt = local_buf_indices.size();
  recvcounts_.resize(atlas::mpi::comm().size());

  atlas::mpi::comm().allGather(sendcnt, recvcounts_.begin(), recvcounts_.end());

  recvdispls_.resize(atlas::mpi::comm().size());
  recvdispls_[0] = 0;
  recvcnt_ = recvcounts_[0];
  for ( size_t jproc = 1; jproc < atlas::mpi::comm().size(); ++jproc ) {
      recvdispls_[jproc] = recvdispls_[jproc - 1] + recvcounts_[jproc - 1];
      recvcnt_ += recvcounts_[jproc];
  }

  ASSERT(recvcnt_ >= static_cast<int>(buffer_indices_->nx()*buffer_indices_->ny()));

  if (myrank == mpiroot) {
    writer_ = std::make_unique<NemoFieldWriter>(file_path, datetimes,
                                                buffer_indices_->nx(), buffer_indices_->ny(),
                                                depths);
    unsorted_buffer_indices_.resize(recvcnt_);
  }

  atlas::mpi::comm().gatherv(&local_buf_indices.front(), local_buf_indices.size(),
                             &unsorted_buffer_indices_.front(),
                             &recvcounts_.front(), &recvdispls_.front(), mpiroot);

  this->write_dimensions();
}
}  // namespace orcamodel
