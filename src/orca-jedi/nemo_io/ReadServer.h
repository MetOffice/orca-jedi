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
#include "atlas/grid.h"  // IWYU pragma: keep
#include "atlas-orca/grid/OrcaGrid.h"  // IWYU pragma: keep

#include "orca-jedi/nemo_io/AtlasIndex.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"
#include "orca-jedi/nemo_io/IoPool.h"
#include "orca-jedi/nemo_io/Redistributor.h"
#include "orca-jedi/nemo_io/IoBackend.h"

#include "eckit/log/Timer.h"
#include "eckit/system/ResourceUsage.h"

namespace orcamodel {
class ReadServer {
 public:
  explicit ReadServer(std::shared_ptr<eckit::Timer> eckit_timer,
      const eckit::PathName& file_path,
      const atlas::Mesh& mesh,
      bool parallel_io = false, size_t n_io_ranks = 0);
  ReadServer(ReadServer &&) = default;
  ReadServer(const ReadServer &) = delete;
  ReadServer &operator=(ReadServer &&) = delete;
  ReadServer &operator=(const ReadServer &) = delete;
  template<class T> void read_var(const std::string& var_name,
      const size_t t_index,
      atlas::array::ArrayView<T, 2>& field_view);
  template<class T> void read_vertical_var(const std::string& var_name,
      atlas::array::ArrayView<T, 2>& field_view);
  size_t get_nearest_datetime_index(const util::DateTime& datetime) const;
  template<class T> T read_fillvalue(const std::string& nemo_var_name) const;

 private:
void log_status() const {
  oops::Log::trace() << "orcamodel::log_status " << eckit_timer_->elapsed() << " "
      << static_cast<double>(eckit::system::ResourceUsage().maxResidentSetSize()) / 1.0e+9
      << " Gb" << std::endl;
}
  template<class T> void read_var_on_root(const std::string& var_name,
      const size_t t_index,
      const size_t z_index,
      std::vector<T>& buffer) const;
  template<class T> void read_vertical_var_on_root(const std::string& var_name,
      const size_t n_levels,
      std::vector<T>& buffer) const;
  template<class T> void fill_field(const std::vector<T>& buffer,
      const size_t z_index,
      atlas::array::ArrayView<T, 2>& field_view) const;
  template<class T> void fill_vertical_field(const std::vector<T>& buffer,
      atlas::array::ArrayView<T, 2>& field_view) const;

  /// \brief Build the parallel I/O pool, redistributor and (on I/O ranks) the
  ///        parallel-netCDF read backend. Only used when parallel_ is true.
  void setup_parallel(const eckit::PathName& file_path, size_t n_io_ranks);
  /// \brief Read one horizontal slice on the pool and scatter it onto the
  ///        per-node field view (fills all nodes, including ghost nodes).
  template<class T> void read_var_parallel(const std::string& var_name,
      const size_t t_index, atlas::array::ArrayView<T, 2>& field_view);

  const size_t mpiroot = 0;
  const size_t myrank = atlas::mpi::rank();
  const atlas::Mesh& mesh_;
  std::unique_ptr<AtlasIndexToBufferIndex> buffer_indices_;
  std::unique_ptr<NemoFieldReader> reader_;
  std::shared_ptr<eckit::Timer> eckit_timer_;

  bool parallel_ = false;
  std::unique_ptr<IoPool> pool_;
  std::unique_ptr<Redistributor> redist_;
  std::unique_ptr<FieldReadBackend> backend_;
};
}  // namespace orcamodel
