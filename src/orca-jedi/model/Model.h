/*
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <ostream>
#include <string>

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
  oops::RequiredParameter<oops::Variables> variables{"model variables", this};
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
    : tstep_(conf.getString("tstep")), geom_(geom), vars_(conf, "model variables") {}

  Model(const Geometry & geom, const Parameters_ & params)
    : tstep_(params.tstep.value()), geom_(geom), vars_(params.variables.value()) {}

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
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const {}
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel
