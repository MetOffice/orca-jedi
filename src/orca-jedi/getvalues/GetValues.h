/* 
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */ 

#pragma once

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation.h"
#include "atlas/functionspace.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"


namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace orcamodel {
  class State;
  class Geometry;

atlas::functionspace::PointCloud atlasObsFuncSpaceFactory (const ufo::Locations & locs);

class GetValues : public util::Printable, private util::ObjectCounter<GetValues> {

 public:
  static const std::string classname() {return "orcamodel::GetValues";}

  GetValues(const Geometry & geom, const ufo::Locations & locs,
            const eckit::Configuration & conf);

  virtual ~GetValues() {};

  void fillGeoVaLs(const State& state, const util::DateTime& dt_begin, 
                   const util::DateTime& dt_end, ufo::GeoVaLs& geovals) const ;

 private:
  void print(std::ostream &) const override {};
  atlas::functionspace::PointCloud atlasObsFuncSpace_;
  atlas::Interpolation interpolator_;
};

}  // namespace orcamodel




