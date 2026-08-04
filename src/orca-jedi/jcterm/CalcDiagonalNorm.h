/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include "eckit/config/Configuration.h"
#include "oops/util/abor1_cpp.h"

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
    ABORT("CalcDiagonalNorm::CalcDiagonalNorm not implemented.");
  }
};

}  // namespace orcamodel

