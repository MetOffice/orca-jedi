/*
 * (C) British Crown Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "orca-jedi/geometry/Geometry.h"

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"

#include "eckit/mpi/Comm.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"


namespace orcamodel {

const std::vector<std::string> Geometry::surface_names(
    {"iiceconc", "sst", "sic_tot_var"});
const std::vector<std::string> Geometry::depth_names({"votemper"});

oops::Variables orcaVariableFactory(const eckit::Configuration & config) {
  eckit::LocalConfiguration nemo_var_mapping;
  config.get("nemo names", nemo_var_mapping);

  std::vector<int> channels{};
  std::vector<std::string> names{};
  for (std::string std_name : nemo_var_mapping.keys()) {
    if (std::find(names.begin(), names.end(), std_name) == names.end()) {
      names.push_back(nemo_var_mapping.getString(std_name));
    }
  }

  return oops::Variables(names, channels);
}

// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm) :
                      comm_(comm), vars_(orcaVariableFactory(config)),
                      variance_vars_(config, "variance names"),
                      grid_(config.getString("grid name")),
                      n_levels_(config.getInt("number levels"))
{
    config.get("nemo names", nemo_var_config);
    auto meshgen_config = grid_.meshgenerator();
    atlas::MeshGenerator meshgen(meshgen_config);
    auto partitioner_config = grid_.partitioner();
    partitioner_ = atlas::grid::Partitioner(partitioner_config);
    mesh_ = meshgen.generate(grid_, partitioner_.partition(grid_));

    funcSpace_ = atlas::functionspace::NodeColumns(
        mesh_, atlas::option::halo(0));
    //    mesh_, atlas::option::halo(this->source_mesh_halo()));
}

// -----------------------------------------------------------------------------
Geometry::~Geometry() {}

// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const
{
  // Number of geoval levels in vars
  std::vector<size_t> varSizes(vars.size());
  std::fill(varSizes.begin(), varSizes.end(), 0);

  for (size_t i=0; i < vars.size(); ++i) {
    if (std::find(surface_names.begin(), surface_names.end(),
          nemo_var_name(vars[i])) != surface_names.end()) {
      varSizes[i] = 1;
    } else if (std::find(depth_names.begin(), depth_names.end(),
          nemo_var_name(vars[i])) != depth_names.end()) {
      varSizes[i] = n_levels_;
    } else {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \" ";
      err_stream << "\" " << vars[i] << " not recognised. " << std::endl;
      throw eckit::NotImplemented(err_stream.str(), Here());
    }
  }
  return varSizes;
}

// -----------------------------------------------------------------------------
atlas::FunctionSpace * Geometry::atlasFunctionSpace() const {
  std::string err_message =
    "orcamodel::Geometry::atlasFunctionSpace Not implemented ";
  throw eckit::NotImplemented(err_message, Here());
  atlas::FunctionSpace* r = nullptr;
  return r;
}

atlas::FieldSet * Geometry::atlasFieldSet() const {
  std::string err_message =
    "orcamodel::Geometry::atlasFieldSet Not implemented ";
  throw eckit::NotImplemented(err_message, Here());
  atlas::FieldSet* r = nullptr;
  return r;
}

const oops::Variables & Geometry::variables() const {
  return vars_;
}

const bool Geometry::variable_in_variable_type(std::string variable_name,
  std::string variable_type) const {
  bool is_bkg_var = variance_vars_.has(variable_name);
  if (variable_type == "background variance") {
    return is_bkg_var;
  } else {
    return !is_bkg_var;
  }
}

void Geometry::print(std::ostream & os) const {
  os << "Not Implemented";
}


}  // namespace orcamodel
