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
      const std::vector<double>& lats,
      const std::vector<double>& lons) {
      size_t nlocs = lats.size();

      // Setup observation functionspace
      oops::Log::trace() << "orcamodel::Interpolator:: creating "
                         << "atlasObsFuncSpace with nlocs = " << nlocs
                         << std::endl;
      atlas::Field points("lonlat", atlas::array::make_datatype<double>(),
          atlas::array::make_shape(nlocs, 2));
      auto arrv_t = atlas::array::make_view<double, 2>(points);
      for (unsigned int j = 0; j < nlocs; ++j) {
        arrv_t(j, 1) = lats[j];
        arrv_t(j, 0) = lons[j];
      }
      oops::Log::trace() << "orcamodel::Interpolator:: creating "
                         << "atlasObsFuncSpace ... done" << std::endl;

      return atlas::functionspace::PointCloud(std::move(points));
  }

  Interpolator::Interpolator(const eckit::Configuration & conf,
      const Geometry & geom, const std::vector<double>& lats,
      const std::vector<double>& lons) :
      nlocs_(lats.size()),
      atlasObsFuncSpace_(atlasObsFuncSpaceFactory(lats, lons)),
      interpolator_(eckit::LocalConfiguration(conf, "atlas-interpolator"),
                    geom.functionSpace(),
                    atlasObsFuncSpace_),
      comm_(geom.getComm()) {
    params_.validateAndDeserialize(conf);
    oops::Log::trace() << "orcamodel::Interpolator:: conf:" << conf
                       << std::endl;
    if (nlocs_ == 0) {
      oops::Log::trace() << "orcamodel::Interpolator:: nlocs == 0" << std::endl;
    }
  }

  void Interpolator::apply(const oops::Variables& vars, const State& state,
      const std::vector<bool> & mask,
      std::vector<double>& result) const {
    oops::Log::trace() << "[" << comm_.rank()
                       << "] orcamodel::Interpolator::apply starting "
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
    oops::Log::debug() << "orcamodel::Interpolator::apply nvars = " << nvars
                       << " nlocs_ = " << nlocs_;
    size_t nvals = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      nvals += nlocs_ * varSizes[jvar];
      oops::Log::debug() << " varSizes[" << jvar << "] = " << varSizes[jvar];
    }
    result.resize(nvals);
    oops::Log::debug() << " nvals = " << nvals << std::endl;

    std::size_t out_idx = 0;
    for (size_t jvar=0; jvar < nvars; ++jvar) {
      auto gv_varname = vars[jvar];
      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>(
          atlas::option::name(gv_varname) |
          atlas::option::levels(varSizes[jvar]));
      interpolator_.execute(state.stateFields()[gv_varname], tgt_field);
      auto field_view = atlas::array::make_view<double, 2>(tgt_field);
      atlas::field::MissingValue mv(state.stateFields()[gv_varname]);
      bool has_mv = static_cast<bool>(mv);
      for (std::size_t klev=0; klev < varSizes[jvar]; ++klev) {
        for (std::size_t iloc=0; iloc < nlocs_; iloc++) {
          if (has_mv && mv(field_view(iloc, klev))) {
            result[out_idx] = util::missingValue(result[out_idx]);
          } else {
            result[out_idx] = field_view(iloc, klev);
          }
          ++out_idx;
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
