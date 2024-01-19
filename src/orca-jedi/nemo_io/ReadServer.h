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

namespace orcamodel {
class ReadServer {
 public:
  explicit ReadServer(const eckit::PathName& file_path,
      const atlas::Mesh& mesh);
  ReadServer(ReadServer &&) = default;
  ReadServer(const ReadServer &) = default;
  ReadServer &operator=(ReadServer &&) = default;
  ReadServer &operator=(const ReadServer &) = default;
  void read_var(const std::string& var_name,
      const size_t t_index,
      atlas::array::ArrayView<double, 2>& field_view);
  void read_vertical_var(const std::string& var_name,
      atlas::array::ArrayView<double, 2>& field_view);
  size_t get_nearest_datetime_index(const util::DateTime& datetime) const;
  template<class T> T read_fillvalue(const std::string& nemo_var_name) const;

 private:
  void read_var_on_root(const std::string& var_name,
      const size_t t_index,
      const size_t z_index,
      std::vector<double>& buffer) const;
  void read_vertical_var_on_root(const std::string& var_name,
      const size_t n_levels,
      std::vector<double>& buffer) const;
  void fill_field(const std::vector<double>& buffer,
      const size_t z_index,
      atlas::array::ArrayView<double, 2>& field_view) const;
  void fill_vertical_field(const std::vector<double>& buffer,
      atlas::array::ArrayView<double, 2>& field_view) const;
  const size_t mpiroot = 0;
  const size_t myrank = atlas::mpi::rank();
  const atlas::Mesh& mesh_;
  const OrcaIndex index_glbarray_;
  std::unique_ptr<NemoFieldReader> reader_;
};
}  // namespace orcamodel
