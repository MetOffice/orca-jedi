/*
 * (C) British Crown Copyright 2020-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation.h"
#include "atlas/functionspace.h"
#include "atlas/field/MissingValue.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

#include "orca-jedi/getvalues/GetValues.h"

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

  atlas::functionspace::PointCloud atlasObsFuncSpaceFactory(
      const ufo::Locations & locs) {
      size_t nlocs = locs.size();
      if (nlocs == 0) {
        std::stringstream err_stream;
        err_stream << "orcamodel::GetValues::"
                   << " Constructor called with no locations " << std::endl;
        err_stream << "         "
                   << "this might mean that there are no observations "
                   << "within the time window" << std::endl;
        err_stream << locs << std::endl;
        throw eckit::BadValue(err_stream.str(), Here());
      }

      // Setup observation functionspace
      std::vector<atlas::PointXY> atlasPoints(nlocs);
      oops::Log::trace() << "orcamodel::GetValues:: atlasPoints with nlocs = "
                         << nlocs << std::endl;
      for (atlas::idx_t i=0; i < nlocs; ++i) {
        atlasPoints[i] = atlas::PointXY(locs.lons()[i], locs.lats()[i]);
      }

      return atlas::functionspace::PointCloud(atlasPoints);
  }

  GetValues::GetValues(const Geometry & geom, const ufo::Locations & locs,
            const eckit::Configuration & conf) :
      atlasObsFuncSpace_(atlasObsFuncSpaceFactory(locs)),
      interpolator_(eckit::LocalConfiguration(conf, "atlas-interpolator"),
                    geom.funcSpace(),
                    atlasObsFuncSpace_ ) {
    oops::Log::trace() << "orcamodel::GetValues:: conf:" << conf
                       << std::endl;
    oops::Log::debug() << "orcamodel::GetValues:: atlasObsFuncSpace_:"
                       << atlasObsFuncSpace_ << std::endl;
  }

  void GetValues::fillGeoVaLs(const State& state,
      const util::DateTime& dt_begin, const util::DateTime& dt_end,
      ufo::GeoVaLs& geovals) const {
    std::size_t nlocs = geovals.nlocs();

    oops::Log::trace() << "orcamodel::GetValues::fillGeoVaLs starting "
                       << std::endl;

    std::vector<double> vals(nlocs);
    size_t nvars = geovals.getVars().size();
    for (size_t j=0; j < nvars; ++j) {
      auto gv_varname = geovals.getVars()[j];
      if (!state.variables().has(gv_varname)) {
        std::stringstream err_stream;
        err_stream << "orcamodel::GetValues::fillGeoVals geovals varname \" "
                   << "\" " << gv_varname << " not found in the model state. "
                   << std::endl;
        err_stream << "    add the variable to the state variables and "
                   << "add a mapping from the geometry to that variable."
                   << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }

      auto nemo_var_name = state.geometry()->nemo_var_name(gv_varname);

      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
          atlas::option::name(gv_varname));
      interpolator_.execute(state.stateFields()[nemo_var_name], tgt_field);
      auto field_view = atlas::array::make_view<double, 1>(tgt_field);
      atlas::field::MissingValue mv(state.stateFields()[nemo_var_name]);
      bool has_mv = static_cast<bool>(mv);
      for (std::size_t i=0; i < nlocs; i++) {
        if (has_mv && mv(field_view(i))) {
          vals[i] = util::missingValue(field_view(i));
        } else {
          vals[i] = field_view(i);
        }
      }

      geovals.putAtLevel(vals, gv_varname, 0);
    }
    oops::Log::debug() << "orcamodel::GetValues::fillGeoVaLs geovals print: "
                       << geovals << std::endl;
    oops::Log::trace() << "orcamodel::GetValues::fillGeoVaLs done "
                       << std::endl;
  }

}  // namespace orcamodel




