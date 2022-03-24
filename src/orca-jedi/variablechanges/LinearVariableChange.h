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
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/variablechanges/LinearVariableChangeParameters.h"

namespace orcamodel {

// -----------------------------------------------------------------------------

class LinearVariableChange : public util::Printable {
 public:
  typedef LinearVariableChangeParameters Parameters_;
  static const std::string classname() {
    return "orcamodel::LinearVariableChange";}

  LinearVariableChange(const Geometry &, const Parameters_ &) {}

  void setTrajectory(const State &, const State &) {}

  void multiply(Increment &, const oops::Variables &) const {}
  void multiplyInverse(Increment &, const oops::Variables &) const {}
  void multiplyAD(Increment &, const oops::Variables &) const {}
  void multiplyInverseAD(Increment &, const oops::Variables &) const {}

 private:
  void print(std::ostream & out) const override {
    out << "orcamodel::LinearVariableChange: Not Implemented"; };
};

}  // namespace orcamodel
