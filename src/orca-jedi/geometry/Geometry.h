/*
 * (C) British Crown Copyright 2017-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once 

#include <map>
#include <iostream>
#include <string>
#include <vector>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"
#include "atlas/runtime/Log.h"

#include "eckit/mpi/Comm.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace atlas {
  class Field;
  class FieldSet;
  class Mesh;

}

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Geometry handles for ORCA model.

  oops::Variables orcaVariableFactory(const eckit::Configuration & config);

class Geometry : public util::Printable {
 public:

  //Geometry(const Parameters_ &, const eckit::mpi::Comm &);
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~Geometry();

  std::vector<size_t> variableSizes(const oops::Variables &) const;
  const eckit::mpi::Comm & getComm() const {return comm_;}

  // SABER interface (stubs at present)
  // Note atlasFunctionSpace is incapable of dealing with horizontal and
  //      vertical stagger.
  atlas::FunctionSpace * atlasFunctionSpace() const;
  atlas::FieldSet * atlasFieldSet() const;

  const oops::Variables & variables() const;

  const atlas::Mesh & mesh() const {return mesh_;};
  const atlas::functionspace::NodeColumns & funcSpace() const {return funcSpace_;};
  const std::string nemo_var_name(std::string std_name) const {return nemo_var_config.getString(std_name);};
  const atlas::idx_t & source_mesh_halo() const {return 0;};
  const bool variable_in_variable_type(std::string variable_name, std::string variable_type) const;

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  oops::Variables vars_;
  oops::Variables variance_vars_;
  eckit::LocalConfiguration nemo_var_config;
  size_t n_levels_;
  atlas::Grid grid_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::functionspace::NodeColumns funcSpace_;

  static const std::vector<std::string> surface_names;
  static const std::vector<std::string> depth_names;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel 
