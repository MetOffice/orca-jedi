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
#include "atlas/functionspace.h"
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

#include "orca-jedi/geometry/GeometryParameters.h"

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
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~Geometry();

  std::vector<size_t> variableSizes(const oops::Variables &) const;
  std::vector<std::string> variableNemoSpaces(const oops::Variables & vars)
      const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  const oops::Variables & variables() const;
  void latlon(std::vector<double> & lats, std::vector<double> & lons,
              const bool halo) const;
  const atlas::FunctionSpace & functionSpace() const {return funcSpace_;}
  const atlas::FieldSet & extraFields() const {return nofields_;}

  const atlas::Grid & grid() const {return grid_;}
  const atlas::Mesh & mesh() const {return mesh_;}
  const std::string nemo_var_name(const std::string std_name) const;
  const bool variable_in_variable_type(std::string variable_name,
    std::string variable_type) const;
  bool levelsAreTopDown() const {return true;}
  std::string distributionType() const {
      return params_.partitioner.value().value_or("serial");}

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  oops::Variables vars_;
  size_t n_levels_;
  OrcaGeometryParameters params_;
  atlas::Grid grid_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::functionspace::NodeColumns funcSpace_;
  atlas::FieldSet nofields_;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel
