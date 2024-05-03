/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/nemo_io/WriteServer.h"

#include<algorithm>

#include "atlas-orca/grid/OrcaGrid.h"


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
  std::vector<T> sorted_buffer(orca_indices_.nx() * orca_indices_.ny());
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
  size_t size = orca_indices_.nx() * orca_indices_.ny();
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
  std::vector<T> buffer;
  size_t size = orca_indices_.nx() * orca_indices_.ny();
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
  atlas::array::ArrayView<double, 2> lonlat{atlas::array::make_view<double, 2>(
                                              mesh_.nodes().lonlat())};
  atlas::field::MissingValue missingValue(mesh_.nodes().lonlat());
  size_t size = orca_indices_.nx() * orca_indices_.ny();

  const std::vector<double> lon_buffer = this->gather_field_on_root(lonlat, 0, missingValue);
  const std::vector<double> lat_buffer = this->gather_field_on_root(lonlat, 1, missingValue);

  if (myrank == mpiroot) {
    ASSERT(lat_buffer.size() == size);
    ASSERT(lon_buffer.size() == size);
    writer_->write_dimensions(lat_buffer, lon_buffer);
  }
}

WriteServer::WriteServer(std::shared_ptr<eckit::Timer> eckit_timer,
    const eckit::PathName& file_path,
    const atlas::Mesh& mesh,
    const std::vector<util::DateTime> datetimes,
    const std::vector<double> depths,
    bool is_serial) : mesh_(mesh),
  orca_indices_(mesh), eckit_timer_(eckit_timer), is_serial_(is_serial),
  n_levels_(depths.size()) {
  oops::Log::trace() << "State(ORCA)::nemo_io::WriteServer::WriteServer" << std::endl;
  std::vector<size_t> local_buf_indices;

  if (is_serial_) {
    if (myrank == mpiroot) {
      writer_ = std::make_unique<NemoFieldWriter>(file_path, datetimes,
                                                  orca_indices_.nx(), orca_indices_.ny(),
                                                  depths);
    }
    this->write_dimensions();
    recvcnt_ = orca_indices_.nx() * orca_indices_.ny();
    return;
  }

  auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));
  auto orcaGrid = atlas::OrcaGrid(mesh_.grid());
  for (atlas::idx_t i_node = 0; i_node < mesh_.nodes().size(); ++i_node) {
    const size_t i = ij(i_node, 0);
    const size_t j = ij(i_node, 1);
    local_buf_indices.emplace_back(orca_indices_(i, j));
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

  ASSERT(recvcnt_ >= static_cast<int>(orca_indices_.nx()*orca_indices_.ny()));

  if (myrank == mpiroot) {
    writer_ = std::make_unique<NemoFieldWriter>(file_path, datetimes,
                                                orca_indices_.nx(), orca_indices_.ny(),
                                                depths);
    unsorted_buffer_indices_.resize(recvcnt_);
  }

  atlas::mpi::comm().gatherv(&local_buf_indices.front(), local_buf_indices.size(),
                             &unsorted_buffer_indices_.front(),
                             &recvcounts_.front(), &recvdispls_.front(), mpiroot);

  this->write_dimensions();
}
}  // namespace orcamodel
