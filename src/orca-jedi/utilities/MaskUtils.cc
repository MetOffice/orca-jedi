/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/utilities/MaskUtils.h"

#include <algorithm>

#include "atlas/array/MakeView.h"
#include "atlas/array/DataType.h"
#include "atlas/field/MissingValue.h"

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

namespace orcamodel {

/// \brief Apply a mask to a FieldSet, setting masked points to missing value.
///
/// For each field in the FieldSet, where mask==0, the field value is set to
/// the field's missing_value. Fields without missing_value metadata are skipped.
///
/// The mask may be:
/// - Surface (nNodes, 1): broadcast to all levels of each field.
/// - Volumetric (nNodes, nLevels): applied per-level. If the field has fewer
///   levels than the mask, only the relevant levels are used.
///
/// \param mask    An int32 mask field (shape nNodes x nLevels).
///               Values: 0=masked, 1=valid.
/// \param fields  The FieldSet to mask in-place.
void applyMaskToFields(
  const atlas::Field & mask,
  atlas::FieldSet & fields) {
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  const atlas::idx_t maskLevels = mask_view.shape(1);

  for (atlas::idx_t i = 0; i < fields.size(); ++i) {
    atlas::Field& field = fields[i];
    atlas::field::MissingValue mv(field);
    if (!mv) {
      oops::Log::debug() << "applyMaskToFields: field '" << field.name()
                         << "' has no missing_value metadata, skipping"
                         << std::endl;
      continue;
    }

    const auto applyMask = [&](auto typeVal) {
      using T = decltype(typeVal);
      auto field_view = atlas::array::make_view<T, 2>(field);
      T fill = field.metadata().get<T>("missing_value");
      const atlas::idx_t nNodes = field_view.shape(0);
      const atlas::idx_t nLevels = field_view.shape(1);
      ASSERT(field_view.shape(0) == mask_view.shape(0));

      if (maskLevels == 1) {
        // Surface mask: broadcast level 0 to all field levels
        for (atlas::idx_t j = 0; j < nNodes; ++j) {
          if (mask_view(j, 0) == 0) {
            for (atlas::idx_t k = 0; k < nLevels; ++k) {
              field_view(j, k) = fill;
            }
          }
        }
      } else {
        // Volumetric mask: apply per-level
        const atlas::idx_t levelsToMask = std::min(nLevels, maskLevels);
        for (atlas::idx_t j = 0; j < nNodes; ++j) {
          for (atlas::idx_t k = 0; k < levelsToMask; ++k) {
            if (mask_view(j, k) == 0) {
              field_view(j, k) = fill;
            }
          }
        }
      }
    };

    if (field.datatype() == atlas::array::DataType::real64()) {
      applyMask(double{});
    } else if (field.datatype() == atlas::array::DataType::real32()) {
      applyMask(float{});
    } else {
      oops::Log::warning() << "applyMaskToFields: field '" << field.name()
                           << "' has unsupported datatype, skipping"
                           << std::endl;
    }
  }
  oops::Log::info() << "applyMaskToFields: applied mask to "
                    << fields.size() << " fields" << std::endl;
}

}  // namespace orcamodel
