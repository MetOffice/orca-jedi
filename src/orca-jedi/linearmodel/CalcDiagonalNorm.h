/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

#include "oops/util/Logger.h"

namespace orcamodel {

class CalcDiagonalNorm {
 public:
  static const std::string classname() {return "orcamodel::CalcDiagonalNorm";}

  // Perform calculation in constructor
  CalcDiagonalNorm(const orcamodel::State &,
                   const orcamodel::Geometry &,
                   orcamodel::Increment &,
                   orcamodel::Increment &,
                   const eckit::LocalConfiguration &);
  ~CalcDiagonalNorm();
};

}  // namespace orcamodel
