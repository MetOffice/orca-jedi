/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/interpolator/Interpolator.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field/MissingValue.h"
#include "atlas/functionspace.h"  // IWYU pragma: keep
#include "atlas/interpolation.h"  // IWYU pragma: keep
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"
#include "oops/util/missingValues.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"

namespace ufo {
class GeoVaLs;
}

namespace orcamodel {

atlas::functionspace::PointCloud atlasObsFuncSpaceFactory(
    const std::vector<double>& lats, const std::vector<double>& lons) {
  const size_t nlocs = lats.size();

  // Setup observation functionspace
  oops::Log::trace() << "orcamodel::Interpolator:: creating "
                     << "atlasObsFuncSpace with nlocs = " << nlocs << std::endl;
  atlas::Field points("lonlat", atlas::array::make_datatype<double>(),
                      atlas::array::make_shape(nlocs, 2));
  auto arrv_t = atlas::array::make_view<double, 2>(points);
  for (unsigned int j = 0; j < nlocs; ++j) {
    arrv_t(j, 1) = lats[j];
    arrv_t(j, 0) = lons[j];
  }
  oops::Log::trace() << "orcamodel::Interpolator:: creating "
                     << "atlasObsFuncSpace ... done" << std::endl;
  return atlas::functionspace::PointCloud(std::move(points));
}

Interpolator::Interpolator(const eckit::Configuration& conf,
                           const Geometry& geom,
                           const std::vector<double>& lats,
                           const std::vector<double>& lons)
    : nlocs_(lats.size()),
      atlasObsFuncSpace_(atlasObsFuncSpaceFactory(lats, lons)),
      interpolator_(eckit::LocalConfiguration(conf, "atlas-interpolator"),
                    geom.functionSpace(), atlasObsFuncSpace_),
      comm_(geom.getComm()) {
  params_.validateAndDeserialize(conf);
  oops::Log::trace() << "orcamodel::Interpolator:: conf:" << conf << std::endl;
  if (nlocs_ == 0) {
    oops::Log::trace() << "orcamodel::Interpolator:: nlocs == 0" << std::endl;
  }

  // Extract the atlas-interpolator configuration
  eckit::LocalConfiguration fwd_conf(conf, "atlas-interpolator");

  // Check if adjoint is enabled
  bool has_adjoint = fwd_conf.has("adjoint") && fwd_conf.getBool("adjoint");

  if (has_adjoint) {
    // Adjoint interpolator: copy configuration but remove non_linear setting
    // (adjoint interpolation only works for linear schemes)
    eckit::LocalConfiguration adjoint_conf(fwd_conf);
    fwd_conf.remove("adjoint");

    // Forward interpolator: use configuration as-is (will have non_linear if present)
    interpolator_ = atlas::Interpolation(fwd_conf,
                                         geom.functionSpace(), atlasObsFuncSpace_);

    if (fwd_conf.has("non_linear")) {
      adjoint_conf.remove("non_linear");
      oops::Log::debug() << "orcamodel::Interpolator: Using asymmetric interpolation - "
                        << "forward with non-linear, adjoint without" << std::endl;
    }

    interpolator_adjoint_ = atlas::Interpolation(adjoint_conf,
                                                  geom.functionSpace(), atlasObsFuncSpace_);
  } else {
    // Adjoint not enabled: single interpolator for forward only
    interpolator_ = atlas::Interpolation(fwd_conf,
                                         geom.functionSpace(), atlasObsFuncSpace_);
  }

  // Store land mask from geometry's extra fields if available
  if (geom.extraFields().has("gmask")) {
    gmask_ = geom.extraFields()["gmask"];
    auto gmask_view = atlas::array::make_view<int, 2>(gmask_);
    oops::Log::debug() << "orcamodel::Interpolator: Stored gmask field from geometry" << std::endl;
    oops::Log::debug() << "  gmask shape: (" << gmask_view.shape(0) << ", "
                       << gmask_view.shape(1) << ")" << std::endl;
  } else if (comm_.rank() == 0) {
    oops::Log::debug() << "orcamodel::Interpolator: gmask field not found in geometry"
                       << " and land masking will not be applied" << std::endl;
  }
}

/// \brief Preprocess the data before performing the interpolation.
/// \param fields the atlas FieldSet collection of fields used in the interpolation.
void Interpolator::preprocess(atlas::FieldSet & fields) {
  fields.haloExchange();
}

/// \brief Preprocess the data before performing the adjoint interpolation.
/// \param fields the atlas FieldSet collection of fields used in the interpolation.
void Interpolator::preprocessAD(atlas::FieldSet & fields) {
  fields.adjointHaloExchange();
  fields.set_dirty();
}

void Interpolator::apply(const oops::Variables& vars, const State& state,
                         const std::vector<bool>& mask,
                         std::vector<double>& result) const {
  oops::Log::trace() << "[" << comm_.rank()
                     << "] orcamodel::Interpolator::apply starting "
                     << std::endl;

  const size_t nvars = vars.size();

  for (size_t j = 0; j < nvars; ++j) {
    if (!state.variables().has(vars[j].name())) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Interpolator::apply varname \" "
                 << "\" " << vars[j].name() << " not found in the model state."
                 << std::endl;
      err_stream << "    add the variable to the state variables and "
                 << "add a mapping from the geometry to that variable."
                 << std::endl;
      throw eckit::BadParameter(err_stream.str(), Here());
    }
  }

  const std::vector<size_t> varSizes = state.geometry()->variableSizes(vars);
  oops::Log::debug() << "orcamodel::Interpolator::apply nvars = " << nvars
                     << " nlocs_ = " << nlocs_;
  size_t nvals = 0;
  for (size_t jvar = 0; jvar < nvars; ++jvar) {
    nvals += nlocs_ * varSizes[jvar];
    oops::Log::debug() << " varSizes[" << jvar << "] = " << varSizes[jvar];
  }
  oops::Log::debug() << " nvals = " << nvals << std::endl;
  result.resize(nvals);

  auto res_iter = result.begin();
  for (size_t jvar = 0; jvar < nvars; ++jvar) {
    const auto execute = [&](auto typeVal) {
      using T = decltype(typeVal);
      executeInterpolation<T>(vars[jvar].name(), varSizes[jvar], state, mask,
                              res_iter);
    };

    ApplyForFieldType(execute,
                      state.geometry()->fieldPrecision(vars[jvar].name()),
                      std::string("orcamodel::Interpolator::apply '") +
                          vars[jvar].name() + "' field type not recognised");
  }
  ASSERT(result.size() == nvals);
  oops::Log::trace() << "orcamodel::Interpolator::apply done " << std::endl;
}

/// \brief Execute the atlas interpolation on a single field.
/// \param gv_varname The GeoVaLs variable name of the source field.
/// \param var_size The number of observation locations in the interpolant.
/// \param state The state object containing the model state.
/// \param mask The locations where the interpolant should be set.
/// \param iter Reference to the interator into the output vector.
template <class T>
void Interpolator::executeInterpolation(
    const std::string& gv_varname, const size_t var_size, const State& state,
    const std::vector<bool>& mask, std::vector<double>::iterator& iter) const {
  atlas::Field tgt_field = atlasObsFuncSpace_.createField<T>(
      atlas::option::name(gv_varname) | atlas::option::levels(var_size));
  interpolator_.execute(state.stateFields()[gv_varname], tgt_field);
  auto field_view = atlas::array::make_view<T, 2>(tgt_field);
  atlas::field::MissingValue mv(state.stateFields()[gv_varname]);
  bool has_mv = static_cast<bool>(mv);
  for (std::size_t iloc = 0; iloc < nlocs_; iloc++) {
    if (mask[iloc]) {
      for (std::size_t klev = 0; klev < var_size; ++klev) {
        if (has_mv && mv(field_view(iloc, klev))) {
          *iter = util::missingValue<double>();
        } else {
          *iter = static_cast<double>(field_view(iloc, klev));
        }
        std::advance(iter, 1);
      }
    } else {
      // skip past elements not required according to the mask.
      std::advance(iter, var_size);
    }
  }
}

template void Interpolator::executeInterpolation<double>(
    const std::string& gv_varname, const size_t var_size, const State& state,
    const std::vector<bool>& mask, std::vector<double>::iterator& iter) const;
template void Interpolator::executeInterpolation<float>(
    const std::string& gv_varname, const size_t var_size, const State& state,
    const std::vector<bool>& mask, std::vector<double>::iterator& iter) const;

/// \brief Interpolate from model space to observation space
/// \param vars Oops variables
/// \param inc Increment object (input)
/// \param mask Mask vector in observation space
/// \param result Result (output) vector in observation space
void Interpolator::apply(const oops::Variables& vars, const Increment& inc,
                         const std::vector<bool>& mask,
                         std::vector<double>& result) const {
  const size_t nvars = vars.size();

  for (size_t j = 0; j < nvars; ++j) {
    if (!inc.variables().has(vars[j])) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Interpolator::apply varname \" "
                 << "\" " << vars[j] << " not found in the model increment."
                 << std::endl;
      err_stream << "    add the variable to the increment variables and "
                 << "add a mapping from the geometry to that variable."
                 << std::endl;
      throw eckit::BadParameter(err_stream.str(), Here());
    }
  }

  const std::vector<size_t> varSizes = inc.geometry()->variableSizes(vars);
  size_t nvals = 0;
  for (size_t jvar = 0; jvar < nvars; ++jvar) nvals += nlocs_ * varSizes[jvar];
  result.resize(nvals);

  std::size_t out_idx = 0;
  for (size_t jvar = 0; jvar < nvars; ++jvar) {
    auto gv_varname = vars[jvar].name();

    // Create source field for interpolation - if gmask available, create a masked copy
    // Otherwise use the increment field directly
    atlas::Field src_field;
    const double mv = util::missingValue<double>();

    if (gmask_) {
      // Create a copy of the increment field with land points set to missing value
      // This way we don't modify the original increment field
      const atlas::Field& inc_field = inc.incrementFields()[gv_varname];
      src_field = inc.geometry()->functionSpace().createField<double>(
          atlas::option::name(gv_varname) |
          atlas::option::levels(varSizes[jvar]));

      auto src_view = atlas::array::make_view<double, 2>(src_field);
      auto inc_view = atlas::array::make_view<double, 2>(inc_field);
      auto gmask_view = atlas::array::make_view<int, 2>(gmask_);

      size_t total_points = 0;
      size_t land_points = 0;
      size_t nonzero_ocean = 0;
      double min_ocean = 1e30;
      double max_ocean = -1e30;

      for (atlas::idx_t jnode = 0; jnode < inc_view.shape(0); ++jnode) {
        for (atlas::idx_t klev = 0; klev < inc_view.shape(1); ++klev) {
          total_points++;
          // gmask is 2D surface field (27118, 1), so always use level 0
          if (gmask_view(jnode, klev) == 0) {
            // Land point: set to missing value
            land_points++;
            src_view(jnode, klev) = mv;
          } else {
            // Ocean point: copy value from increment
            src_view(jnode, klev) = inc_view(jnode, klev);
            if (std::abs(inc_view(jnode, klev)) > 1e-12) {
              nonzero_ocean++;
              min_ocean = std::min(min_ocean, inc_view(jnode, klev));
              max_ocean = std::max(max_ocean, inc_view(jnode, klev));
            }
          }
        }
      }

      oops::Log::debug() << "orcamodel::Interpolator::apply(increment): Variable " << gv_varname
                         << std::endl;
      oops::Log::debug() << "  Total points: " << total_points
                         << ", Land points: " << land_points
                         << " (" << (100.0 * land_points / total_points) << "%)" << std::endl;
      oops::Log::debug() << "  Non-zero ocean points: " << nonzero_ocean << std::endl;
      if (nonzero_ocean > 0) {
        oops::Log::debug() << "  Ocean value range: [" << min_ocean << ", " << max_ocean << "]"
                           << std::endl;
      }
    } else {
      // No mask available, use increment field directly
      src_field = inc.incrementFields()[gv_varname];
      if (comm_.rank() == 0) {
        oops::Log::debug() << "orcamodel::Interpolator::apply(increment): No gmask available for "
                           << gv_varname << std::endl;
      }
    }

    atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
        atlas::option::name(gv_varname) |
        atlas::option::levels(varSizes[jvar]));
    interpolator_.execute(src_field, tgt_field);
    auto field_view = atlas::array::make_view<double, 2>(tgt_field);
    atlas::field::MissingValue field_mv(src_field);
    bool has_mv = static_cast<bool>(field_mv);
    for (std::size_t iloc = 0; iloc < nlocs_; iloc++) {
      for (std::size_t klev = 0; klev < varSizes[jvar]; ++klev) {
        if (has_mv && field_mv(field_view(iloc, klev))) {
          result[out_idx] = util::missingValue<double>();
        } else {
          result[out_idx] = field_view(iloc, klev);
        }
        ++out_idx;
      }
    }
  }
}

/// \brief Interpolate from observation space to model space
/// \param vars Oops variables
/// \param inc Increment object (output)
/// \param mask Mask (observation space) vector
/// \param resultin Values (observation space) vector (input)
void Interpolator::applyAD(const oops::Variables& vars, Increment& inc,
                           const std::vector<bool>& mask,
                           const std::vector<double>& resultin) const {
  oops::Log::trace() << "orcamodel::Interpolator::applyAD start " << std::endl;

  const size_t nvars = vars.size();

  for (size_t j = 0; j < nvars; ++j) {
    if (!inc.variables().has(vars[j])) {
      std::stringstream err_stream;
      err_stream << "orcamodel::Interpolator::apply varname \" "
                 << "\" " << vars[j] << " not found in the model increment."
                 << std::endl;
      err_stream << "    add the variable to the increment variables and "
                 << "add a mapping from the geometry to that variable."
                 << std::endl;
      throw eckit::BadParameter(err_stream.str(), Here());
    }
  }

  const std::vector<size_t> varSizes = inc.geometry()->variableSizes(vars);

  std::size_t out_idx = 0;
  for (size_t jvar = 0; jvar < nvars; ++jvar) {
    auto gv_varname = vars[jvar].name();
    auto tgt_field = atlasObsFuncSpace_.createField<double>(
        atlas::option::name(gv_varname) |
        atlas::option::levels(varSizes[jvar]));

    // Use the OOPS missing value as the default for increment fields,
    // however in the future we might need to switch to a value defined in
    // configuration or read from file, therefore we keep it logically separate
    // for now.
    const double inc_default_missing_value = util::missingValue<double>();
    tgt_field.metadata().set("missing_value", inc_default_missing_value);
    tgt_field.metadata().set("missing_value_type", "approximately-equals");
    tgt_field.metadata().set("missing_value_epsilon", 1e-6);

    // Copying observation array vector to an atlas observation field
    // (tgt_field)
    auto field_view = atlas::array::make_view<double, 2>(tgt_field);

    for (std::size_t iloc = 0; iloc < nlocs_; iloc++) {
      for (std::size_t klev = 0; klev < varSizes[jvar]; ++klev) {
        // Only apply adjoint for valid (non-masked) observations
        // For masked or missing observations, set to ZERO (not missing value)
        // because:
        // 1. Zero = no gradient contribution (additive identity)
        // 2. The increment field already has proper missing value structure
        // 3. execute_adjoint doesn't respect missing value metadata, so we need zero
        if (!mask[iloc]) {
          // Masked observation: set to zero (no contribution to adjoint)
          field_view(iloc, klev) = 0.0;
        } else if (resultin[out_idx] == util::missingValue<double>() ||
                   !std::isfinite(resultin[out_idx])) {
          // Missing/invalid value in observation: set to zero (no contribution to adjoint)
          field_view(iloc, klev) = 0.0;
        } else {
          // Valid observation: use the value
          field_view(iloc, klev) = resultin[out_idx];
        }
        ++out_idx;
      }
    }

    // halo exchange update ghost points
    std::shared_ptr<const Geometry> geom = inc.geometry();
    geom->functionSpace().haloExchange(inc.incrementFields()[gv_varname]);

    // BEFORE adjoint: Apply land mask to increment field if gmask is available
    // Set land points to 0.0 to prevent adjoint from spreading gradients over land
    if (gmask_) {
      atlas::Field& inc_field = inc.incrementFields()[gv_varname];
      auto inc_view = atlas::array::make_view<double, 2>(inc_field);
      auto gmask_view = atlas::array::make_view<int, 2>(gmask_);

      size_t total_points = 0;
      size_t land_points = 0;
      size_t zeroed_before = 0;
      double min_before = 1e30;
      double max_before = -1e30;

      for (atlas::idx_t jnode = 0; jnode < inc_view.shape(0); ++jnode) {
        for (atlas::idx_t klev = 0; klev < inc_view.shape(1); ++klev) {
          total_points++;
          // gmask == 0 means land (masked), gmask == 1 means ocean (valid)
          if (gmask_view(jnode, klev) == 0) {
            land_points++;
            if (std::abs(inc_view(jnode, klev)) > 1e-12) {
              zeroed_before++;
            }
            inc_view(jnode, klev) = 0.0;
          } else {
            min_before = std::min(min_before, inc_view(jnode, klev));
            max_before = std::max(max_before, inc_view(jnode, klev));
          }
        }
      }
      oops::Log::debug() << "orcamodel::Interpolator::applyAD: Variable " << gv_varname
                         << " BEFORE adjoint" << std::endl;
      oops::Log::debug() << "  Total points: " << total_points
                         << ", Land points: " << land_points
                         << " (" << (100.0 * land_points / total_points) << "%)" << std::endl;
      oops::Log::debug() << "  Zeroed non-zero land points: " << zeroed_before << std::endl;
      oops::Log::debug() << "  Ocean value range before adjoint: [" << min_before << ", "
                         << max_before << "]" << std::endl;
    }

    // Use the adjoint-specific interpolator (without non-linear settings)
    interpolator_adjoint_.execute_adjoint(inc.incrementFields()[gv_varname], tgt_field);

    // AFTER adjoint: Apply land mask again to zero out any spurious values
    if (gmask_) {
      atlas::Field& inc_field = inc.incrementFields()[gv_varname];
      auto inc_view = atlas::array::make_view<double, 2>(inc_field);
      auto gmask_view = atlas::array::make_view<int, 2>(gmask_);

      size_t zeroed_after = 0;
      double min_ocean_after = 1e30;
      double max_ocean_after = -1e30;
      size_t nonzero_ocean = 0;

      for (atlas::idx_t jnode = 0; jnode < inc_view.shape(0); ++jnode) {
        for (atlas::idx_t klev = 0; klev < inc_view.shape(1); ++klev) {
          if (gmask_view(jnode, klev) == 0) {
            // Land point
            if (std::abs(inc_view(jnode, klev)) > 1e-12) {
              zeroed_after++;
            }
            inc_view(jnode, klev) = 0.0;
          } else {
            // Ocean point
            if (std::abs(inc_view(jnode, klev)) > 1e-12) {
              nonzero_ocean++;
              min_ocean_after = std::min(min_ocean_after, inc_view(jnode, klev));
              max_ocean_after = std::max(max_ocean_after, inc_view(jnode, klev));
            }
          }
        }
      }
      oops::Log::debug() << "orcamodel::Interpolator::applyAD: AFTER adjoint" << std::endl;
      oops::Log::debug() << "  Zeroed non-zero land points: " << zeroed_after << std::endl;
      oops::Log::debug() << "  Non-zero ocean points: " << nonzero_ocean << std::endl;
      if (nonzero_ocean > 0) {
        oops::Log::debug() << "  Ocean value range after adjoint: ["
                         << min_ocean_after << ", " << max_ocean_after << "]" << std::endl;
      }
    }
  }  // jvar

  oops::Log::trace() << "orcamodel::Interpolator::applyAD done " << std::endl;
}

void Interpolator::print(std::ostream& os) const {
  os << "orcamodel::Interpolator: " << std::endl;
  os << "  Obs function space " << atlasObsFuncSpace_ << std::endl;
  os << "  Interpolator " << interpolator_ << std::endl;
}

}  // namespace orcamodel
