/*
 * (C) British Crown Copyright 2023 Met Office
 */

#pragma once

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace orcamodel {
  class ModelBias;
  class ModelBiasIncrement;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasCovariance : public util::Printable,
                       private eckit::NonCopyable,
                       private util::ObjectCounter<ModelBiasCovariance> {
 public:
  static const std::string classname()
                           {return "orcamodel::ModelBiasCovariance";}

/// Constructor, destructor
  ModelBiasCovariance(const eckit::Configuration & conf,
                      const Geometry &): conf_(conf) {}
  ~ModelBiasCovariance() {}

/// Linear algebra operators
  void linearize(const ModelBias &, const Geometry &) {}
  void multiply(const ModelBiasIncrement &,
                ModelBiasIncrement) const {}
  void inverseMultiply(const ModelBiasIncrement &,
                       ModelBiasIncrement) const {}
  void randomize(ModelBiasIncrement &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel

