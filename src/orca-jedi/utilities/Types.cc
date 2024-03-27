/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/utilities/Types.h"

namespace orcamodel {

template <typename T> const netCDF::NcType NetCDFTypeMap<T>::ncType = netCDF::ncDouble;
template <> const netCDF::NcType NetCDFTypeMap<float>::ncType = netCDF::ncFloat;

template const netCDF::NcType NetCDFTypeMap<float>::ncType;
template const netCDF::NcType NetCDFTypeMap<double>::ncType;

}  // namespace orcamodel
