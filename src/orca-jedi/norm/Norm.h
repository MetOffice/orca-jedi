/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"

namespace orcamodel {

class Norm : public util::Printable,
             private util::ObjectCounter<Norm> {
 public:
  static const std::string classname() {return "orcamodel::Norm";}

// Constructor, destructor
  Norm(const Geometry &,
       const oops::Variables &,
       const util::DateTime &,
       const eckit::Configuration &);
  ~Norm();

// Compute values for norm
  void calculate(const State &);

// Apply norm to increment
  void apply(Increment &) const;
  void applyInverse(Increment &) const;

 private:
  void print(std::ostream &) const override;
};

}  // namespace orcamodel
