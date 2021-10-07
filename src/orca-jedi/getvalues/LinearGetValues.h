/*
 * (C) British Crown Copyright 2017-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"


#include "oops/base/VariableChangeBase.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"

// ----------------------------------------------------------------------------

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace orcamodel {
  class Increment;
  class State;
  class Geometry;
  struct UnifiedModelTraits;

// -----------------------------------------------------------------------------

class LinearGetValues : public util::Printable,
  private util::ObjectCounter<LinearGetValues> {
 public:
  static const std::string classname() {return "orcamodel::LinearGetValues";}

  LinearGetValues(const Geometry &, const ufo::Locations &,
      const eckit::Configuration &);
  virtual ~LinearGetValues();

  void setTrajectory(const State & state, const util::DateTime & t1,
      const util::DateTime & t2, ufo::GeoVaLs & geovals) {}
  void fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1,
      const util::DateTime & t2, ufo::GeoVaLs & geovals) const {}
  void fillGeoVaLsAD(Increment & inc, const util::DateTime & t1,
      const util::DateTime & t2, const ufo::GeoVaLs & geovals) const {}

 private:
  void print(std::ostream &) const;
  ufo::Locations locs_;
  std::shared_ptr<const Geometry> geom_;
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel
