/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/utilities/OrcaModelTraits.h"
#include "oops/base/VariableChangeBase.h"

namespace orcamodel {

// -----------------------------------------------------------------------------

class VarChaModel2GeoVaLs: public oops::VariableChangeBase<OrcaModelTraits>,
                           private util::ObjectCounter<VarChaModel2GeoVaLs> {
 public:
  static const std::string classname() {
    return "orcamodel::VarChaModel2GeoVaLs";
  }
  VarChaModel2GeoVaLs(const Geometry &, const eckit::Configuration &);
  ~VarChaModel2GeoVaLs();
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  std::shared_ptr<const Geometry> geom_;
  void print(std::ostream &) const override;
};

}  // namespace orcamodel
