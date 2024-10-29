/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <tuple>

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/utilities/Types.h"

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/parallel/mpi/mpi.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "eckit/mpi/Comm.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/system/ResourceUsage.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace {
/// \brief Construct an atlas grid given a string containing either a grid name
///        or a path to a grid specification yaml configuration file.
/// \param[in]     grid_specification  string containing the path/name.
/// \return        constructed atlas grid object.
atlas::Grid construct_grid_from_name(std::string grid_specification) {
  auto grid_name = grid_specification;
  eckit::PathName grid_spec_path(grid_specification);

  std::vector<std::string> orca_grid_names;
  for (auto && orca_type : std::vector<std::string>{"F", "T", "U", "V", "W"}) {
    for (auto && resolution : std::vector<std::string>{"1", "2", "025", "12"}) {
      orca_grid_names.emplace_back("ORCA" + resolution + "_" + orca_type);
      orca_grid_names.emplace_back("eORCA" + resolution + "_" + orca_type);
    }
  }
  auto grid = atlas::Grid();
  if (std::find(std::begin(orca_grid_names), std::end(orca_grid_names), grid_name)
      != std::end(orca_grid_names)) {
    grid = atlas::Grid{grid_name};
  } else if (grid_spec_path.exists()) {
    grid = atlas::Grid{atlas::Grid::Spec{grid_spec_path}};
  } else {
    std::stringstream err_stream;
    err_stream << "orcamodel::Geometry:: grid  \"" << grid_specification
               << "\" " << " is neither a valid named grid,"
               << " nor a path to a grid specification. " << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  return grid;
}
}  // namespace

namespace orcamodel {

oops::Variables orcaVariableFactory(const eckit::Configuration & config) {
  OrcaGeometryParameters params;
  params.validateAndDeserialize(config);

  oops::Variables variables{};
  std::vector<std::string> names{};
  for (const NemoFieldParameters& nemoVariable :
        params.nemoFields.value()) {
    std::string name = nemoVariable.name.value();
    if (std::find(names.begin(), names.end(), name) == names.end()) {
      names.emplace_back(name);
      variables.push_back(oops::Variable(name));
    }
  }

  return variables;
}

// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & config,
                   const eckit::mpi::Comm & comm) :
                      comm_(comm), vars_(orcaVariableFactory(config)),
                      n_levels_(config.getInt("number levels")),
                      eckit_timer_(new eckit::Timer("Geometry(ORCA): ", oops::Log::trace()))
{
    eckit_timer_->start();
    log_status();
    params_.validateAndDeserialize(config);

    grid_ = construct_grid_from_name(params_.gridName.value());

    int64_t halo = params_.sourceMeshHalo.value();
    std::string partitioner_name = params_.partitioner.value();
    if ( ( (partitioner_name == "serial") || (comm.size() == 1) )
         && (halo > 0) ) {
      halo = 0;
      partitioner_name = "serial";
      oops::Log::info() << "Warning: forcing halo = 0 and serial partitioner"
                        << " as settings imply all processors have all data" << std::endl;
    }
    auto meshgen_config = grid_.meshgenerator()
                          | atlas::option::halo(halo);

    atlas::MeshGenerator meshgen(meshgen_config);
    log_status();
    auto partitioner_config = grid_.partitioner();
    partitioner_config.set("type", partitioner_name);
    partitioner_ = atlas::grid::Partitioner(partitioner_config);
    log_status();
    mesh_ = meshgen.generate(grid_, partitioner_);
    log_status();
    funcSpace_ = atlas::functionspace::NodeColumns(
        mesh_, atlas::option::halo(halo));
    log_status();

    if (params_.extraFieldsInit.value().value_or(false)) {
      // Fill extra geometry fields for BUMP / SABER
      // these are area, vunit, hmask, gmask
      extraFields_ = atlas::FieldSet();

      oops::Log::debug() << "calling orcagrid" << std::endl;
      atlas::OrcaGrid orcaGrid = mesh_.grid();
      nx_ = orcaGrid.nx() + orcaGrid.haloWest() + orcaGrid.haloEast();
      ny_ = orcaGrid.ny() + orcaGrid.haloNorth();

      oops::Log::debug() << "orcagrid nx_ " << nx_ << " ny_ " << ny_
                         << std::endl;

      // Vertical unit - DJL set to something sensible?
      atlas::Field vunit = funcSpace_.createField<double>(
        atlas::option::name("vunit") | atlas::option::levels(n_levels_));
      auto field_view = atlas::array::make_view<double, 2>(vunit);
      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
           field_view(j, k) = 1.;
        }
      }
      extraFields_->add(vunit);

      // halo mask / owned
      atlas::Field hmask = funcSpace_.createField<int32_t>(
        atlas::option::name("owned") | atlas::option::levels(n_levels_));
      auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());

      auto field_view1 = atlas::array::make_view<int32_t, 2>(hmask);
      for (atlas::idx_t j = 0; j < field_view1.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view1.shape(1); ++k) {
          int x, y;
          std::tie(x, y) = xypt(j, nx_);
          // 0 mask, 1 ocean
  //        if (ghost(j) || x >= 181 || y >= 147 ) {field_view1(j, k) = 0;    // DJL
          if (ghost(j) || x >= nx_-1 || y >= ny_-1 ) {field_view1(j, k) = 0;    // DJL
          } else {field_view1(j, k) = 1;}
        }
      }
      // Add field
      oops::Log::debug() << "Geometry adding hmask (set to all ocean except halo)"
                         << std::endl;
      extraFields_->add(hmask);

      // geometry mask
      atlas::Field gmask = funcSpace_.createField<int32_t>(
        atlas::option::name("gmask") | atlas::option::levels(n_levels_));

      auto field_view2 = atlas::array::make_view<int32_t, 2>(gmask);
      for (atlas::idx_t j = 0; j < field_view2.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view2.shape(1); ++k) {
          int x, y;
          std::tie(x, y) = xypt(j, nx_);
          // DJL hardwired to orca2 needs generalising
  //        if (ghost(j) || x >= 181 || y >= 147 ) {field_view2(j, k) = 0;   // DJL
          if (ghost(j) || x >= nx_-1 || y >= ny_-1 ) {field_view2(j, k) = 0;    // DJL
          // 0 mask, 1 ocean
          } else {field_view2(j, k) = 1;}
        }
      }
      // Add field
      oops::Log::debug() << "Geometry adding gmask (set to all ocean except halo)"
                         << std::endl;
      extraFields_->add(gmask);

      atlas::Field area = funcSpace_.createField<double>(
        atlas::option::name("area") | atlas::option::levels(n_levels_));
      auto field_view3 = atlas::array::make_view<double, 2>(area);

      for (atlas::idx_t j = 0; j < field_view3.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view3.shape(1); ++k) {
          field_view3(j, k) = 4e10;           // DJL should change   2 degrees
        }
      }
      log_status();

      // Add field
      extraFields_->add(area);
    }
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
      if (nemoField.name.value() == vars[i].name()) {
        if (nemoField.modelSpace.value() == "surface") {
          varSizes[i] = 1;
        } else {
          varSizes[i] = n_levels_;
        }
      }
    }
    if (varSizes[i] == 0) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \" ";
      err_stream << "\" " << vars[i].name() << " not recognised. " << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }
  }
  return varSizes;
}

void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
    const bool halo) const {
  const auto lonlat = atlas::array::make_view<double, 2>(funcSpace_.lonlat());
  const auto ghosts = atlas::array::make_view<int32_t, 1>(
      mesh_.nodes().ghost());
  const auto haloDistance = atlas::array::make_view<int32_t, 1>(
      mesh_.nodes().halo());
  auto isRequired = [&](const size_t nodeElem) {
    if (halo) {
      return !ghosts(nodeElem) || (haloDistance(nodeElem) > 0);
    }
    return !ghosts(nodeElem);
  };
  const size_t npts = funcSpace_.size();
  for (size_t nodeElem = 0; nodeElem < npts; ++nodeElem) {
    if (isRequired(nodeElem)) {
      lons.emplace_back(lonlat(nodeElem, 0));
      lats.emplace_back(lonlat(nodeElem, 1));
    }
  }
}

// -----------------------------------------------------------------------------
/// \brief Give the space of nemo field for each variable - surface, volume or
///         vertical. at the moment we need this distinction to read 3D depth
///         data from a 1D array
/// \param[in]     vars  variables to check.
/// \return        vector of variable Nemo model spaces.
std::vector<std::string> Geometry::variableNemoSpaces(
    const oops::Variables & vars) const
{
  std::vector<std::string> varNemoSpaces(vars.size(), "");

  auto nemoFields = params_.nemoFields.value();

  for (size_t i=0; i < vars.size(); ++i) {
    for (const auto & nemoField : nemoFields) {
      if (nemoField.name.value() == vars[i].name()) {
        if (nemoField.modelSpace.value() == "surface" ||
            nemoField.modelSpace.value() == "volume" ||
            nemoField.modelSpace.value() == "vertical" ) {
          varNemoSpaces[i] = nemoField.modelSpace.value();
        } else {
            std::stringstream err_stream;
            err_stream << "orcamodel::Geometry::variableNemoSpaces modelSpace"
                       << " \"" << nemoField.modelSpace.value()
                       << "\" not recognised for field \""
                       << nemoField.name.value() << "\"." << std::endl;
            throw eckit::BadValue(err_stream.str(), Here());
        }
      }
    }
    if (varNemoSpaces[i] == "") {
      std::stringstream err_stream;
      err_stream << "orcamodel::Geometry::variableSizes variable name \"";
      err_stream << vars[i] << "\" not available in the state. ";
      err_stream << "Either add this state variable to the model ";
      err_stream << "configuration or remove the corresponding obs variable";
      err_stream << " from the filters configuration." << std::endl;
      throw eckit::BadValue(err_stream.str(), Here());
    }
  }
  return varNemoSpaces;
}

const oops::Variables & Geometry::variables() const {
  return vars_;
}

/// \brief Check if a variable's data is a member of a type (e.g if it can be
///        sourced from the background file, error file, or MDT file).
/// \param[in]     variable_name  Name of variable.
/// \param[in]     variable_type  Type of variable.
/// \return        Boolean for membership.
const bool Geometry::variable_in_variable_type(std::string variable_name,
  std::string variable_type) const {
  auto nemoFields = params_.nemoFields.value();
  for (const auto & nemoField : nemoFields) {
    if (nemoField.name.value() == variable_name) {
      std::string type = nemoField.variableType.value();
      return type == variable_type;
    }
  }

  std::stringstream err_stream;
  err_stream << "orcamodel::Geometry::variable_in_variable_type variable name ";
  err_stream << "\"" << variable_name << "\" not recognised. " << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

/// \brief Data type of the atlas field holding the named variable data.
/// \param[in]     variable_name  Name of variable.
/// \return        orcamodel::FieldDType enum.
FieldDType Geometry::fieldPrecision(std::string variable_name) const {
  auto nemoFields = params_.nemoFields.value();
  for (const auto & nemoField : nemoFields) {
    if (nemoField.name.value() == variable_name) {
       return nemoField.fieldPrecision.value();
    }
  }

  std::stringstream err_stream;
  err_stream << "orcamodel::Geometry::fieldPrecision variable name ";
  err_stream << "\"" << variable_name << "\" not recognised. " << std::endl;
  throw eckit::BadValue(err_stream.str(), Here());
}

void Geometry::print(std::ostream & os) const {
  os << "Not Implemented";
}

void Geometry::log_status() const {
  oops::Log::trace() << "orcamodel::log_status " << eckit_timer_->elapsed() << " "
      << static_cast<double>(eckit::system::ResourceUsage().maxResidentSetSize()) / 1.0e+9
      << " Gb" << std::endl;
}

void Geometry::set_gmask(atlas::Field & field) const {
  oops::Log::debug() << "Geometry setting gmask from field missing values" << std::endl;
  oops::Log::debug() << "Fieldname: " << field.name() << std::endl;  // DJL

  atlas::Field gmask = extraFields_.field("gmask");

  atlas::field::MissingValue mv(field);
  bool has_mv = static_cast<bool>(mv);
  oops::Log::debug() << "has_mv " << has_mv << std::endl;       // DJL

  const auto setGmaskField = [&](auto typeVal) {
    using T = decltype(typeVal);
    auto field_viewin = atlas::array::make_view<T, 2>(field);
    auto field_viewgm = atlas::array::make_view<int32_t, 2>(gmask);
//    auto lonlat_view = atlas::array::make_view<double, 2>(funcSpace_.lonlat());

    if (has_mv) {
      for (atlas::idx_t j = 0; j < field_viewgm.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_viewgm.shape(1); ++k) {
// only change values that are currently unmasked ( 0 mask, 1 ocean )
          if (field_viewgm(j, k) == 1) {
            if (mv(field_viewin(j, k))) {
              field_viewgm(j, k) = 0;
            }
          }
        }
      }
    }
  };
  ApplyForFieldType(setGmaskField,
                    fieldPrecision(field.name()),
                    std::string("Geometry(ORCA)::set_gmask ")
                    + field.name() + "' field type not recognised");
  log_status();
}

// Determine x,y location from jpt
std::tuple<int, int> xypt(int jpt, int nx) {
  int y = jpt / nx;
  int x = jpt - y*nx;
  return std::make_tuple(x, y); }
}  // namespace orcamodel
