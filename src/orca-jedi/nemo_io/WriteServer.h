/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>
#include <memory>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"  // IWYU pragma: keep
#include "atlas/mesh.h"  // IWYU pragma: keep
#include "atlas/field/MissingValue.h"
#include "atlas-orca/grid/OrcaGrid.h"  // IWYU pragma: keep

#include "orca-jedi/nemo_io/AtlasIndex.h"
#include "orca-jedi/nemo_io/NemoFieldWriter.h"
#include "orca-jedi/nemo_io/IoPool.h"
#include "orca-jedi/nemo_io/Redistributor.h"
#include "orca-jedi/nemo_io/IoBackend.h"

#include "eckit/log/Timer.h"
#include "eckit/system/ResourceUsage.h"  // IWYU pragma: keep

namespace orcamodel {
class WriteServer {
 public:
  explicit WriteServer(std::shared_ptr<eckit::Timer> eckit_timer,
      const eckit::PathName& file_path,
      const atlas::Mesh& mesh,
      const std::vector<util::DateTime> datetimes,
      const std::vector<double> depths, bool is_serial,
      bool parallel_io = false, size_t n_io_ranks = 0);
  WriteServer(WriteServer &&) = default;
  WriteServer(const WriteServer &) = delete;
  WriteServer &operator=(WriteServer &&) = delete;
  WriteServer &operator=(const WriteServer &) = delete;

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

  /// \brief Build the parallel I/O pool, redistributor and (on I/O ranks) the
  ///        parallel-netCDF backend. Only used when parallel_ is true.
  void setup_parallel(const eckit::PathName& file_path,
      const std::vector<util::DateTime>& datetimes,
      const std::vector<double>& depths, size_t n_io_ranks);
  /// \brief Gather a per-node level into this rank's I/O slab (fill masked).
  template<class T> std::vector<T> field_to_slab(
      const atlas::array::ArrayView<T, 2>& field_view, const size_t i_level,
      const atlas::field::MissingValue& missingValue) const;
  void write_dimensions_parallel();
  template<class T> void write_vol_var_parallel(const std::string& var_name,
      const size_t t_index, const atlas::field::MissingValue& missingValue,
      const atlas::array::ArrayView<T, 2>& field_view);
  template<class T> void write_surf_var_parallel(const std::string& var_name,
      const size_t t_index, const atlas::field::MissingValue& missingValue,
      const atlas::array::ArrayView<T, 2>& field_view);

  const size_t mpiroot = 0;
  const size_t myrank = atlas::mpi::rank();
  const atlas::Mesh& mesh_;
  std::unique_ptr<AtlasIndexToBufferIndex> buffer_indices_;
  std::vector<size_t> unsorted_buffer_indices_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;
  std::unique_ptr<NemoFieldWriter> writer_;
  std::shared_ptr<eckit::Timer> eckit_timer_;
  bool is_serial_;
  const size_t n_levels_;
  int recvcnt_;

  bool parallel_ = false;
  std::unique_ptr<IoPool> pool_;
  std::unique_ptr<Redistributor> redist_;
  std::unique_ptr<FieldWriteBackend> backend_;
};
}  // namespace orcamodel
