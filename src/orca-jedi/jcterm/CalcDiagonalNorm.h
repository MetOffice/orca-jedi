/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"

namespace orcamodel {

class CalcDiagonalNorm {
 public:
  CalcDiagonalNorm(const State &,
                   const Geometry &,
                   Increment &,
                   Increment &,
                   const eckit::Configuration &) {
    throw eckit::NotImplemented(Here());
  }
};

}  // namespace orcamodel

