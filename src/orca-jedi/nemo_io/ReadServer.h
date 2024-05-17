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
#include "atlas/grid.h"
#include "atlas-orca/grid/OrcaGrid.h"

#include "orca-jedi/nemo_io/OrcaIndex.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/log/Timer.h"
#include "eckit/system/ResourceUsage.h"

namespace orcamodel {
class ReadServer {
 public:
  static const std::string classname() {return "orcamodel::ReadServer";}
  explicit ReadServer(std::shared_ptr<eckit::Timer> eckit_timer,
      const eckit::PathName& file_path,
      const atlas::Mesh& mesh);
  ReadServer(ReadServer &&) = default;
  ReadServer(const ReadServer &) = default;
  ReadServer &operator=(ReadServer &&) = default;
  ReadServer &operator=(const ReadServer &) = default;
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
  const size_t mpiroot = 0;
  const size_t myrank = atlas::mpi::rank();
  const atlas::Mesh& mesh_;
  const OrcaIndexToBufferIndex orca_buffer_indices_;
  std::unique_ptr<NemoFieldReader> reader_;
  std::shared_ptr<eckit::Timer> eckit_timer_;
};
}  // namespace orcamodel
