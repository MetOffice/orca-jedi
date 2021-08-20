/*
 * (C) British Crown Copyright 2017-2018 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <ostream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/log/CodeLocation.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/Locations.h"
#include "ufo/GeoVaLs.h"

#include "orca-jedi/errorcovariance/ErrorCovariance.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/increment/Increment.h"

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const oops::Variables & vars,
                     const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(), time_(time),
    incrementFields_()
{
  oops::Log::trace() << "Increment(orca)::create and zero 1" << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const Increment & other)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
//  unifiedmodel_increment_change_resol_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "Increment(orca)::constr constructed from other."
                     << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
    oops::Log::trace() << "Increment(orca)::constr create, zero"
                       << "and copy = " << copy
                       << std::endl;
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  oops::Log::trace() << "Increment(orca)::constr create and copy 2"
                     << std::endl;
}
// -----------------------------------------------------------------------------
Increment::~Increment() {
  oops::Log::trace() << "Increment(orca)::destructed" << std::endl;
}

// SABER interface (stubs at present)
// -----------------------------------------------------------------------------
void Increment::setAtlas(atlas::FieldSet * fs) const {
  oops::Log::trace() << "Increment(orca)::setAtlas" << std::endl;
}

//------------------------------------------------------------------------
void Increment::toAtlas(atlas::FieldSet * fs) const {
  oops::Log::trace() << "Increment(orca)::toAtlas" << std::endl;
}

// -----------------------------------------------------------------------------
void Increment::fromAtlas(atlas::FieldSet * fs)  {
  oops::Log::trace() << "Increment(orca)::fromAtlas" << std::endl;
}


}  // namespace orcamodel
