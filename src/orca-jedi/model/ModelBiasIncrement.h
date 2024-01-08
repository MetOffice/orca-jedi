/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <iostream>
#include <vector>

#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

namespace eckit {
  class Configuration;
}

namespace orcamodel {
  class ModelBias;
  class ModelBiasCovariance;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasIncrement : public util::Printable,
                           public util::Serializable{
 public:
/// Constructor, destructor
  ModelBiasIncrement(const Geometry &,
                     const eckit::Configuration &) {}
  ModelBiasIncrement(const ModelBiasIncrement &,
                     const bool) {}
  ModelBiasIncrement(const ModelBiasIncrement &,
                     const eckit::Configuration &) {}
  ~ModelBiasIncrement() {}

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &) {}
  void zero() {}
  ModelBiasIncrement & operator=(const ModelBiasIncrement &) {
    return *this;}
  ModelBiasIncrement & operator+=(const ModelBiasIncrement &) {
    return *this;}
  ModelBiasIncrement & operator-=(const ModelBiasIncrement &) {
    return *this;}
  ModelBiasIncrement & operator*=(const double) {return *this;}
  void axpy(const double, const ModelBiasIncrement &) {}
  double dot_product_with(const ModelBiasIncrement &) const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

/// Serialize and deserialize
  std::size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, std::size_t &) override {}

 private:
  explicit ModelBiasIncrement(const ModelBiasCovariance &);
  void print(std::ostream & os) const override {}
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel

