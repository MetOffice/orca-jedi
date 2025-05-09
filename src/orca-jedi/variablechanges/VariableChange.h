/*
 * (C) British Crown Copyright 2024 Met Office
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
  void changeVar(State & state, const oops::Variables & variables) const {
    state.subsetFieldSet(variables);
  }

  void changeVarInverse(State &, const oops::Variables &) const {}

 private:
  void print(std::ostream & out) const override {
    out << "orcamodel::VariableChange: Not Implemented"; };
};

}  // namespace orcamodel


