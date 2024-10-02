/*
 * (C) Crown Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/Printable.h"

#include "orca-jedi/geometry/Geometry.h"


namespace orcamodel {

class Geometry;


class ModelData : public util::Printable {
 public:
  explicit ModelData(const Geometry &) {}
  ~ModelData() {}

  static const std::string classname() {return "orcamodel::ModelData";}

  const eckit::LocalConfiguration modelData() const;

 private:
  void print(std::ostream &) const override;
};


}  // namespace orcamodel
