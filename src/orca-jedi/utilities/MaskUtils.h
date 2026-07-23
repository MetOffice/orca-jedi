/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include "atlas/field.h"

namespace orcamodel {

/// \brief Apply a mask field to a FieldSet, setting masked points to missing.
///
/// The mask field should be int32 with values 0=masked, 1=valid.
/// - If the mask has 1 level (surface), it is broadcast to all field levels.
/// - If the mask has the same number of levels as the field, per-level masking
///   is applied.
/// Fields without missing_value metadata are skipped.
void applyMaskToFields(
  const atlas::Field & mask,
  atlas::FieldSet & fields);

}  // namespace orcamodel
