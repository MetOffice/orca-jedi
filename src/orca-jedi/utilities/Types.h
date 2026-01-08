/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <netcdf>

#include <string>

#include "eckit/exception/Exceptions.h"

#include "atlas/array/DataType.h"  // IWYU pragma: keep

namespace orcamodel {

//// \brief Enum type for obs variable data types
enum class FieldDType {
    Float,
    Double
};

/// \brief Apply a function for a given FieldDType
template<typename Functor>
void ApplyForFieldType(const Functor& functor, FieldDType field_type,
     const std::string& error_message) {
  if (field_type == FieldDType::Float) {
    functor(float{});
  } else if (field_type == FieldDType::Double) {
    functor(double{});
  } else {
    throw eckit::BadParameter(error_message);
  }
}

/// \brief Apply a function for a given Atlas DataType
template<typename Functor>
void ApplyForFieldType(const Functor& functor, atlas::DataType datatype,
     const std::string& error_message) {
  if (datatype.str() == "real32") {
    functor(float{});
  } else if (datatype.str() == "real64") {
    functor(double{});
  } else {
    throw eckit::BadParameter(error_message);
  }
}

// create a mapping between C++ types and NetCDF type objects
template <typename T>
struct NetCDFTypeMap;

template <>
struct NetCDFTypeMap<float> {
  inline static const netCDF::NcType ncType = netCDF::ncFloat;
};

template <>
struct NetCDFTypeMap<double> {
  inline static const netCDF::NcType ncType = netCDF::ncDouble;
};

template <>
struct NetCDFTypeMap<int> {
  inline static const netCDF::NcType ncType = netCDF::ncInt;
};

}  // namespace orcamodel
