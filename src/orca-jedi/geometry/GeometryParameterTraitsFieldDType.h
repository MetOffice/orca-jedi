/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ParameterTraits.h"

#include "orca-jedi/utilities/Types.h"

namespace orcamodel {

/// Helps with the conversion of FieldDType values to/from strings.
struct FieldDTypeParameterTraitsHelper {
  typedef FieldDType EnumType;
  static constexpr char enumTypeName[] = "FieldDType";
  static constexpr util::NamedEnumerator<EnumType> namedValues[] = {
    { EnumType::Float, "float" },
    { EnumType::Double, "double" }
  };
};

}  // namespace orcamodel

namespace oops {

/// Specialization of ParameterTraits for FieldDType.
template <>
struct ParameterTraits<orcamodel::FieldDType> :
    public EnumParameterTraits<orcamodel::FieldDTypeParameterTraitsHelper>
{};

}  // namespace oops
