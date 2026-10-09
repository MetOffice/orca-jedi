/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "orca-jedi/norm/Norm.h"

#include "eckit/exception/Exceptions.h"

namespace orcamodel {

Norm::Norm(const oops::Variables &,
           const eckit::Configuration &) {
  throw eckit::NotImplemented(Here());
}

void Norm::calculate(const State &) {
  throw eckit::NotImplemented(Here());
}

void Norm::apply(Increment &) const {
  throw eckit::NotImplemented(Here());
}

void Norm::applyInverse(Increment &) const {
  throw eckit::NotImplemented(Here());
}

void Norm::print(std::ostream &) const {
  throw eckit::NotImplemented(Here());
}

}  // namespace orcamodel
