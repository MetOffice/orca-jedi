/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include<exception>

namespace orcamodel {

//// \brief Enum type for obs variable data types
enum class FieldDType {
    Float,
    Double
};

/// \brief Apply a function for a given FieldDType
template<typename Functor>
void ApplyForFieldType(const Functor& functor, FieldDType field_type, std::exception on_fail) {
  if (field_type == FieldDType::Float) {
    functor(float{});
  } else if (field_type == FieldDType::Double) {
    functor(double{});
  } else {
    throw on_fail;
  }
}

}  // namespace orcamodel
