/*
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>
#include <sstream>

#include "eckit/exception/Exceptions.h"
#include "oops/interface/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/DateTime.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/OptionalParameter.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/StateParameters.h"
#include "orca-jedi/utilities/OrcaModelTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace orcamodel {
  class Increment;
  class ModelBias;
  class State;


class OrcaModelParameters : public oops::ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(OrcaModelParameters, ModelParametersBase)

 public:
  /// Model time step
  oops::RequiredParameter<util::Duration> tstep{"tstep", this};
  /// Model variables
  oops::RequiredParameter<std::vector<OrcaStateParameters>> states{
    "states",
    "List of configuration options used to initialize the" +
    " model state at each time step",
    this};
};


// -----------------------------------------------------------------------------
/// OrcaModel model definition.
/*!
 *  OrcaModel nonlinear model definition and configuration parameters.
 */

class Model: public oops::interface::ModelBase<OrcaModelTraits>,
             private util::ObjectCounter<Model> {
 public:
  typedef OrcaModelParameters Parameters_;

  static const std::string classname() {return "orcamodel::Model";}

  Model(const Geometry & geom, const eckit::Configuration & conf)
    : tstep_(conf.getString("tstep")), geom_(geom) {
    parameters_.validateAndDeserialize(conf);
    checkTimeStep();
  }

  Model(const Geometry & geom, const Parameters_ & params)
    : parameters_(params), tstep_(params.tstep.value()), geom_(geom)
       {
    oops::Log::trace() << classname() << "constructor begin" << std::endl;
    checkTimeStep();
    oops::Log::trace() << classname() << "constructor end" << std::endl;
  }

  ~Model() {}

/// Prepare model integration
  void initialize(State &) const {}

/// Model integration
  void step(State &, const ModelBias &) const {}
  int saveTrajectory(State &, const ModelBias &) const { return 0;}

/// Finish model integration
  void finalize(State &) const {}

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}

 private:
  void checkTimeStep() {
    auto stateParameters = parameters_.states.value();
    for (size_t iState = 0; iState < stateParameters.size() - 1; ++iState) {
      const util::DateTime time = stateParameters[iState].date.value();
      const util::DateTime timeNext = stateParameters[iState + 1].date.value();
      std::cout << classname() << "::checkTimeStep state difference "
                << (timeNext - time) << " tstep " << tstep_
                << std::endl;
      if ((timeNext - time).toSeconds() != tstep_.toSeconds()) {
        std::ostringstream msg;
        msg << "state difference " << (timeNext - time)
            << " != tstep " << tstep_ << std::endl;
        throw eckit::UserError(msg.str(), Here());
      }
    }
  }
  void print(std::ostream &) const {}
  util::Duration tstep_;
  Parameters_ parameters_;
  const Geometry geom_;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel
