/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#include <string>

#include "eckit/exception/Exceptions.h"

#include "orca-jedi/variablechanges/LinearVariableChange.h"

#include "oops/util/Printable.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/variablechanges/LinearVariableChangeParameters.h"

namespace orcamodel {

// -----------------------------------------------------------------------------
LinearVariableChange::LinearVariableChange(const Geometry &, const Parameters_ &) {}
LinearVariableChange::LinearVariableChange(const Geometry &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
// LinearVariableChange::~LinearVariableChange() {}

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars)
const {}

void LinearVariableChange::changeVarInverseTL(Increment & dx, const oops::Variables & vars)
const {}

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars)
const {}

void LinearVariableChange::changeVarInverseAD(Increment & dx, const oops::Variables & vars)
const {}

void LinearVariableChange::changeVarTraj(const State & xbg, const oops::Variables & xfg)
{}

}  // namespace orcamodel
