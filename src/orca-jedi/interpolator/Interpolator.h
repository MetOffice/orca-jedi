/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation.h"
#include "atlas/functionspace.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/interpolator/InterpolatorParameters.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"


namespace eckit {
class Configuration;
}

namespace orcamodel {
class State;
class Geometry;
class Increment;

atlas::functionspace::PointCloud atlasObsFuncSpaceFactory(
    const std::vector<double> & locs);

class Interpolator : public util::Printable,
  private util::ObjectCounter<Interpolator> {
 public:
  static const std::string classname() {return "orcamodel::Interpolator";}

  // Parameters will not work properly until support added to oops getvalues
  // typedef OrcaInterpolatorParameters Parameters_;

  Interpolator(const eckit::Configuration & conf, const Geometry & geom,
      const std::vector<double>& lats, const std::vector<double>& lons);

  virtual ~Interpolator() {}

  void apply(const oops::Variables& vars, const State& state,
             const std::vector<bool> & mask,
             std::vector<double>& result) const;
  void apply(const oops::Variables& vars, const Increment& inc,
             const std::vector<bool> & mask,
             std::vector<double>& result) const {
    throw eckit::NotImplemented("Increment interpolation not implemented",
                                Here());
  }
  void applyAD(const oops::Variables& vars, const Increment& inc,
               const std::vector<bool> & mask,
               const std::vector<double> &) const {
    throw eckit::NotImplemented("Adjoint interpolation not implemented",
                                Here());
  }

 private:
  template<class T> void executeInterpolation(
      const std::string& gv_varname,
      size_t var_size,
      const State& state,
      const std::vector<bool> & mask,
      std::vector<double>::iterator& result) const;
  void print(std::ostream &) const override;
  int64_t nlocs_;
  atlas::functionspace::PointCloud atlasObsFuncSpace_;
  atlas::Interpolation interpolator_;
  // Parameters_ params_;
  OrcaInterpolatorParameters params_;
  const eckit::mpi::Comm & comm_;
};

}  // namespace orcamodel
