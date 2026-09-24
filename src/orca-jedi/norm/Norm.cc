/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "orca-jedi/norm/Norm.h"

#include "eckit/exception/Exceptions.h"

namespace orcamodel {

Norm::Norm(const Geometry & resol,
           const oops::Variables & vars,
           const util::DateTime & validTime,
           const eckit::Configuration & conf) {
  throw eckit::NotImplemented(Here());
}

Norm::~Norm() {}

void Norm::calculate(const State & xx) {
  throw eckit::NotImplemented(Here());
}

void Norm::apply(Increment & dx) const {
  throw eckit::NotImplemented(Here());
}

void Norm::applyInverse(Increment & dx) const {
  throw eckit::NotImplemented(Here());
}

void Norm::print(std::ostream & os) const {
  throw eckit::NotImplemented(Here());
}

}  // namespace orcamodel
