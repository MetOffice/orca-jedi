/*
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/interface/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

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

// -----------------------------------------------------------------------------
/// OrcaModel model definition.
/*!
 *  OrcaModel nonlinear model definition and configuration parameters.
 */

class Model: public oops::interface::ModelBase<OrcaModelTraits>,
             private util::ObjectCounter<Model> {
 public:
  static const std::string classname() {return "orcamodel::Model";}

  Model(const Geometry & geom, const eckit::Configuration & conf)
    : tstep_(0), geom_(geom), vars_(conf, "model variables") {}
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
