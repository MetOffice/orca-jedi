/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/variablechanges/VariableChangeParameters.h"

// Forward declarations
namespace oops {
  class Variables;
}

namespace orcamodel {

// -----------------------------------------------------------------------------

class VariableChange: public util::Printable,
  private util::ObjectCounter<VariableChange> {
 public:
  static const std::string classname() {
    return "orcamodel::VariableChange";
  }
  VariableChange(const VariableChangeParameters &, const Geometry &) {}
  VariableChange(const eckit::Configuration &, const Geometry &) {}
  void changeVar(State &, const oops::Variables &) const {}
  void changeVarInverse(State &, const oops::Variables &) const {}

 private:
  void print(std::ostream & out) const override {
    out << "orcamodel::VariableChange: Not Implemented"; };
};

}  // namespace orcamodel


