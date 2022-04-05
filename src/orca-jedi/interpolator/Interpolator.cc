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
#include <utility>

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

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

#include "orca-jedi/interpolator/Interpolator.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
}

namespace orcamodel {
  class State;
  class Geometry;

  atlas::functionspace::PointCloud atlasObsFuncSpaceFactory(
      const std::vector<double>& locs) {
      size_t nlocs = locs.size() / 2;
      if (nlocs == 0) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::"
                   << " Constructor called with no locations " << std::endl;
        err_stream << "         "
                   << "this might mean that there are no observations "
                   << "within the time window" << std::endl;
        err_stream << locs << std::endl;
        throw eckit::BadValue(err_stream.str(), Here());
      }

      // Setup observation functionspace
      oops::Log::trace() << "orcamodel::Interpolator:: creating atlasObsFuncSpace "
                         << "with nlocs = " << nlocs << std::endl;
      atlas::Field points("lonlat", atlas::array::make_datatype<double>(),
          atlas::array::make_shape(nlocs, 2));
      auto arrv_t = atlas::array::make_view<double, 2>(points);
      for ( unsigned int j = 0, latIndex=0, lonIndex=1;
        j < nlocs;
        ++j, latIndex+=2, lonIndex+=2 ) {
        arrv_t(j, 1) = locs[latIndex];
        arrv_t(j, 0) = locs[lonIndex];
      }
      oops::Log::trace() << "orcamodel::Interpolator:: creating atlasObsFuncSpace "
                         << "... done" << std::endl;

      return atlas::functionspace::PointCloud(std::move(points));
  }

  Interpolator::Interpolator(const eckit::Configuration & conf, const Geometry & geom,
      const std::vector<double>& locs) :
      nlocs_(locs.size() / 2),
      atlasObsFuncSpace_(std::move(atlasObsFuncSpaceFactory(locs))),
      interpolator_(eckit::LocalConfiguration(conf, "atlas-interpolator"),
                    geom.funcSpace(),
                    atlasObsFuncSpace_ ) {
    params_.validateAndDeserialize(conf);
    oops::Log::trace() << "orcamodel::Interpolator:: conf:" << conf
                       << std::endl;
    oops::Log::debug() << "orcamodel::Interpolator:: atlasObsFuncSpace_:"
                       << atlasObsFuncSpace_ << std::endl;
  }

  void Interpolator::apply(const oops::Variables& vars, const State& state,
      std::vector<double>& result) const {

    oops::Log::trace() << "orcamodel::Interpolator::apply starting "
                       << std::endl;

    const size_t nvars = vars.size();
    for (size_t j=0; j < nvars; ++j) {
      if (!state.variables().has(vars[j])) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::apply varname \" "
                   << "\" " << vars[j]
                   << " not found in the model state." << std::endl;
        err_stream << "    add the variable to the state variables and "
                   << "add a mapping from the geometry to that variable."
                   << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }
    }
    const std::vector<size_t> varSizes =
      state.geometry()->variableSizes(vars);
    for (size_t j=0; j < nvars; ++j) {
      auto gv_varname = vars[j];
      if (varSizes[j] != 1) {
        std::stringstream err_stream;
        err_stream << "orcamodel::Interpolator::apply interpolating "
                   << "data with levels > 1 not implemented." << std::endl;
        throw eckit::NotImplemented(err_stream.str(), Here());
      }
      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
          atlas::option::name(gv_varname));
      interpolator_.execute(state.stateFields()[gv_varname], tgt_field);
      auto field_view = atlas::array::make_view<double, 1>(tgt_field);
      atlas::field::MissingValue mv(state.stateFields()[gv_varname]);
      bool has_mv = static_cast<bool>(mv);
      for (std::size_t i=0; i < nlocs_; i++) {
        if (has_mv && mv(field_view(i))) {
          result[i] = util::missingValue(field_view(i));
        } else {
          result[i] = field_view(i);
        }
      }
    }
    oops::Log::trace() << "orcamodel::Interpolator::apply done "
                       << std::endl;
  }

  void Interpolator::print(std::ostream & os) const {
    os << "orcamodel::Interpolator: " << std::endl;
    os << "  Obs function space " << atlasObsFuncSpace_ << std::endl;
    os << "  Interpolator " << interpolator_ << std::endl;
  }

}  // namespace orcamodel
