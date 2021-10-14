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
#include "oops/util/Logger.h"


namespace orcamodel {

oops::Variables orcaVariableFactory(const eckit::Configuration & config) {
  Geometry::Parameters__ params;
  params.validateAndDeserialize(config);

  std::vector<int> channels{};
  std::vector<std::string> names{};
  for (const NemoFieldParameters& nemoVariable :
        params.nemoFields.value()) {
    std::string  name = nemoVariable.name.value();
    if (std::find(names.begin(), names.end(), name) == names.end()) {
      names.push_back(name);
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
    params_.validateAndDeserialize(config);
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

const std::string Geometry::nemo_var_name(const std::string std_name) const {
  for (const auto & nemoField : params_.nemoFields.value()) {
    if (std_name == nemoField.name.value()) return nemoField.nemoName.value();
  }
  std::stringstream err_stream;
  err_stream << "orcamodel::Geometry::nemo_var_name variable name \" ";
  err_stream << "\" " << std_name << " not recognised. " << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

// -----------------------------------------------------------------------------
/// \brief Give the number of levels for each provided level - surface variables
///        have 1 level, volumetric variables have "number levels" levels.
/// \param[in]     vars  variables to check. 
/// \return        vector of number of levels in each variable. 
std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const
{
  std::vector<size_t> varSizes(vars.size());
  std::fill(varSizes.begin(), varSizes.end(), 0);

  auto nemoFields = params_.nemoFields.value();

  for (size_t i=0; i < vars.size(); ++i) {
    for (const auto & nemoField : nemoFields) {
      if (nemoField.name.value() == vars[i]) {
        if (nemoField.type.value() == "surface") {
          varSizes[i] = 1;
        } else {
          varSizes[i] = n_levels_;
        }
      }
    }
    if (varSizes[i] == 0) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \" ";
      err_stream << "\" " << vars[i] << " not recognised. " << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
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

/// \brief Check if a variable's data is a member of a type (e.g if it can be
///        sourced from the background file, variance file, or MDT file).
/// \param[in]     variable_name  Name of variable. 
/// \param[in]     variable_type  Type of variable. 
/// \return        Boolean for membership. 
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
