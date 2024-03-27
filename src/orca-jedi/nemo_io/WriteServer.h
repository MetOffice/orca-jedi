/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>
#include <vector>
#include <memory>

#include "oops/util/DateTime.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"
#include "atlas/mesh.h"
#include "atlas/field/MissingValue.h"
#include "atlas-orca/grid/OrcaGrid.h"

#include "orca-jedi/nemo_io/OrcaIndex.h"
#include "orca-jedi/nemo_io/NemoFieldWriter.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/log/Timer.h"
#include "eckit/system/ResourceUsage.h"

namespace orcamodel {
class WriteServer {
 public:
  explicit WriteServer(std::shared_ptr<eckit::Timer> eckit_timer,
      const eckit::PathName& file_path,
      const atlas::Mesh& mesh,
      const std::vector<util::DateTime> datetimes,
      const std::vector<double> depths, bool is_serial);
  WriteServer(WriteServer &&) = default;
  WriteServer(const WriteServer &) = default;
  WriteServer &operator=(WriteServer &&) = default;
  WriteServer &operator=(const WriteServer &) = default;

  template<class T> void write_vol_var(const std::string& var_name,
      const size_t t_index,
      const atlas::field::MissingValue& missingValue,
      const atlas::array::ArrayView<T, 2>& field_view);
  template<class T> void write_surf_var(const std::string& var_name,
      const size_t t_index,
      const atlas::field::MissingValue& missingValue,
      const atlas::array::ArrayView<T, 2>& field_view);

 private:
  void log_status() const {
    oops::Log::trace() << "orcamodel::log_status " << eckit_timer_->elapsed() << " "
        << static_cast<double>(eckit::system::ResourceUsage().maxResidentSetSize()) / 1.0e+9
        << " Gb" << std::endl;
  }
  void write_dimensions();
  template<class T> void write_vol_var_on_root(const std::string& var_name,
      const size_t t_index, const std::vector<T>& buffer);
  template<class T> void write_surf_var_on_root(const std::string& var_name,
      const size_t t_index, const std::vector<T>& buffer);
  template<class T> std::vector<T> gather_field_on_root(
      const atlas::array::ArrayView<T, 2>& field_view, const size_t i_level,
      const atlas::field::MissingValue& missingValue) const;
  template<class T> std::vector<T> sort_buffer(
    const std::vector<T> & buffer) const;
  const size_t mpiroot = 0;
  const size_t myrank = atlas::mpi::rank();
  const atlas::Mesh& mesh_;
  const OrcaIndexToBufferIndex orca_indices_;
  std::vector<size_t> unsorted_buffer_indices_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;
  std::unique_ptr<NemoFieldWriter> writer_;
  std::shared_ptr<eckit::Timer> eckit_timer_;
  bool is_serial_;
  const size_t n_levels_;
  int recvcnt_;
};
}  // namespace orcamodel
