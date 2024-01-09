/*
 * (C) British Crown Copyright 2024 Met Office
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
  std::string err_message =
      "orcamodel::Increment::constructor not implemented";
  throw eckit::NotImplemented(err_message, Here());
}
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const Increment & other)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::string err_message =
      "orcamodel::Increment::constructor(geom, other) not implemented";
  throw eckit::NotImplemented(err_message, Here());
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::string err_message =
      "orcamodel::Increment::constructor(other, copy) not implemented";
  throw eckit::NotImplemented(err_message, Here());
}
// -----------------------------------------------------------------------------
Increment::Increment(const Increment & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::string err_message =
      "orcamodel::Increment::copy constructor not implemented";
  throw eckit::NotImplemented(err_message, Here());
}


}  // namespace orcamodel
