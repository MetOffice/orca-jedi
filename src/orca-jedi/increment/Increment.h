/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <ostream>
#include <string>
#include <memory>
#include <vector>

#include "atlas/field/FieldSet.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/dot_product.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

#include "eckit/exception/Exceptions.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
}

namespace oops {
  class Variables;
  class UnstructuredGrid;
}

namespace orcamodel {
  class Geometry;
  class ModelBiasIncrement;
  class ErrorCovariance;
  class State;

/// orcaModel Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. Similarly visa-versa.
 *  The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

// -----------------------------------------------------------------------------

class Increment : public util::Printable,
                  public util::Serializable,
                  private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "orcamodel::Increment";}

/// Constructor, destructor
  Increment(const Geometry &,
            const oops::Variables &,
            const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  Increment(const Increment &);

/// Basic operators
  void diff(const State &, const State &) {}
  void zero() {}
  void zero(const util::DateTime &) {}
  void ones() {}
  // Increment & operator =(const Increment &);
  // Increment & operator+=(const Increment &);
  // Increment & operator-=(const Increment &);
  // Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true) {}
  double dot_product_with(const Increment &) const {return 0.0; }
  void schur_product_with(const Increment &) {}
  void random() {}
  void dirac(const eckit::Configuration &) {}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0; }
  void print(std::ostream & os) const override {os << "Not Implemented";}

  void updateTime(const util::Duration & dt) {time_ += dt;}

  std::vector<double> rmsByLevel(const std::string &) const {
    ABORT("rmsByLevel not implemented");
    return {};
  }

/// Serialize and deserialize
  std::size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, std::size_t &) override {}

/// Other
  void accumul(const double &, const State &) {}


/// Utilities
  std::shared_ptr<const Geometry> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  const atlas::FieldSet & incrementFields() const {return incrementFields_;}
  atlas::FieldSet & incrementFields() {return incrementFields_;}

  const oops::Variables & variables() const {return vars_;}


/// Data
 private:
  // void print(std::ostream &) const override;
  std::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  util::DateTime time_;
  atlas::FieldSet incrementFields_;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel

