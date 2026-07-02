/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/regridder/Regridder.h"

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

namespace orcamodel {

Regridder::Regridder(const eckit::Configuration& conf,
                     const atlas::FunctionSpace& sourceFunctionSpace,
                     const atlas::FunctionSpace& targetFunctionSpace)
    : sourceFunctionSpace_(sourceFunctionSpace),
      targetFunctionSpace_(targetFunctionSpace),
      interpolation_(conf, sourceFunctionSpace_, targetFunctionSpace_) {}

atlas::FieldSet Regridder::execute(const atlas::FieldSet& source) const {
  atlas::FieldSet target;
  for (atlas::idx_t i = 0; i < source.size(); ++i) {
    target.add(execute(source[i]));
  }
  return target;
}

atlas::Field Regridder::execute(const atlas::Field& source) const {
  atlas::Field target = targetFunctionSpace_.createField(
      atlas::option::name(source.name()) |
      atlas::option::levels(source.levels()) |
      atlas::option::datatype(source.datatype()));
  interpolation_.execute(source, target);
  return target;
}

}  // namespace orcamodel
