/*
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "atlas/trans/Trans.h"
#include "atlas/trans/ifs/TransIFS.h"

#include "eckit/config/Configuration.h"
#include "eckit/memory/NonCopyable.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "orca-jedi/geometry/Geometry.h"

// Forward declarations
namespace atlas {
class FieldSet;
}

namespace oops {
class Variables;
}

namespace orcamodel {
class Increment;
class State;

/// Background error covariance matrix
class ErrorCovariance : public util::Printable,
        private eckit::NonCopyable,
        private util::ObjectCounter<ErrorCovariance> {
 public:
  static const std::string classname() {return "orcamodel::ErrorCovariance";}

  ErrorCovariance(const Geometry &, const oops::Variables &,
    const eckit::Configuration &, const State &, const State &);
  ~ErrorCovariance();

  void linearize(const State &, const Geometry &);
  void multiply(const Increment &, Increment &) const;
  void inverseMultiply(const Increment &, Increment &) const;
  void randomize(Increment &) const;

 private:
  void print(std::ostream &) const;
  std::shared_ptr<const Geometry> geom_;

  util::DateTime time_;
};

}  // namespace orcamodel
