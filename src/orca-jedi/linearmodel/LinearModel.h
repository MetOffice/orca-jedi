/*
 * (C) British Crown Copyright 2025 Met Office
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
class Configuration;
}

namespace util {
class Duration;
}

namespace orcamodel {
class Geometry;
class Increment;
class ModelBias;
class ModelBiasIncrement;
class State;

// -----------------------------------------------------------------------------

/// OrcaModel linearmodel definition.
/*
 *  Empty OrcaModel linear model definition to satisfy oops API.
 */

class LinearModel: public util::Printable,
                   private util::ObjectCounter<LinearModel> {
 public:
  static const std::string classname() {return "orcamodel::LinearModel";}
  static std::vector<std::string> names() {return {"empty-orcaLinearModel"};}

  LinearModel(const Geometry &, const eckit::Configuration &) {
    ABORT("LinearModel not implemented");
  }
  ~LinearModel() {}

/// Model trajectory computation
  void setTrajectory(const State &,
                     State &,
                     const ModelBias &) {}

/// Run TLM and its adjoint
  void initializeTL(Increment &) const  {}
  void stepTL(Increment &, const ModelBiasIncrement &) const  {}
  void finalizeTL(Increment &) const  {}

  void initializeAD(Increment &) const  {}
  void stepAD(Increment &, ModelBiasIncrement &) const  {}
  void finalizeAD(Increment &) const  {}

/// Other utilities
  const util::Duration & timeResolution() const {return emptyDuration_;}
  const util::Duration & stepTrajectory() const {return emptyDuration_;}

 private:
  void print(std::ostream &) const {}

// Data
  const util::Duration emptyDuration_;
};
// -----------------------------------------------------------------------------
}  // namespace orcamodel
