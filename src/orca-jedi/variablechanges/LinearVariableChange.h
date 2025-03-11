/*
 * (C) British Crown Copyright 2024 Met Office
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
#include "orca-jedi/variablechanges/LinearVariableChangeParameters.h"

namespace orcamodel {

// -----------------------------------------------------------------------------

class LinearVariableChange: public util::Printable {
 public:
  typedef LinearVariableChangeParameters Parameters_;
  static const std::string classname() {
    return "orcamodel::LinearVariableChange";
  }

  LinearVariableChange(const Geometry &, const Parameters_ &);
  LinearVariableChange(const Geometry &, const eckit::Configuration &);

  void changeVarTL(Increment &, const oops::Variables &) const;
  void changeVarInverseTL(Increment &, const oops::Variables &) const;
  void changeVarAD(Increment &, const oops::Variables &) const;
  void changeVarInverseAD(Increment &, const oops::Variables &) const;

  void changeVarTraj(const State &, const oops::Variables &);


 private:
  void print(std::ostream & out) const override {
    out << "orcamodel::LinearVariableChange::print Not Implemented"; };
};

}  // namespace orcamodel


