/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "eckit/config/LocalConfiguration.h"

#include "nemovar_jedi/GeometryNV.h"
#include "nemovar_jedi/VariablesNV.h"
#include "nemovar_jedi/ErrorCovarianceNV.h"

#include "orca-jedi/geometry/Geometry.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

#include "orca-jedi/saber_covariance/SaberOrcaJediFilter.h"
// include "saber/blocks/SaberBlockParametersBase.h"
// include "saber/blocks/SaberOuterBlockBase.h"
#include "saber/blocks/SaberCentralBlockBase.h"

namespace saber {
namespace orcajedifilter {

// -----------------------------------------------------------------------------

class OrcaJediFilterParameters : public SaberBlockParametersBase {
  OOPS_CONCRETE_PARAMETERS(OrcaJediFilterParameters, SaberBlockParametersBase)

 public:
    /// Whether to normalize as a localization function
    // oops::Parameter<bool> normalize{"normalize filter variance", true, this};
    /// Filter specifications
    oops::Parameter<eckit::LocalConfiguration> function{"function",
                                                        eckit::LocalConfiguration(), this};
    oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
    oops::RequiredParameter<std::string> covtype{"covtype", this};

    oops::Variables mandatoryActiveVars() const override {return oops::Variables();}
};

// -----------------------------------------------------------------------------

class OrcaJediFilter : public SaberCentralBlockBase {
 public:
  static const std::string classname() {return "saber::orcajedifilter::OrcaJediFilter";}

  typedef OrcaJediFilterParameters Parameters_;

  OrcaJediFilter(const oops::GeometryData &,
                           const oops::Variables &,
                           const eckit::Configuration &,
                           const Parameters_ &,
                           const oops::FieldSet3D &,
                           const oops::FieldSet3D &);

  virtual ~OrcaJediFilter();
//  virtual ~OrcaJediFilter() = default;
//
//  const oops::GeometryData & innerGeometryData() const override {return innerGeometryData_;}
//  const oops::Variables & innerVars() const override {return innerVars_;}

  void multiply(oops::FieldSet3D &) const override;
  void randomize(oops::FieldSet3D &) const override;

  size_t ctlVecSize() const override {return ctlVecSize_;}

//  oops::FieldSet3D generateInnerFieldSet(const oops::GeometryData & innerGeometryData,
//                                        const oops::Variables & innerVars) const override;

//  oops::FieldSet3D generateOuterFieldSet(const oops::GeometryData & outerGeometryData,
//                                        const oops::Variables & outerVars) const override;

//  std::shared_ptr<const nv::GeometryNV> getNVgeometryPtr() const {return nvgeom_;}

 private:
  void print(std::ostream &) const override;

  /// Parameters
  Parameters_ params_;
  /// Active variables
  const oops::Variables activeVars_;
  const oops::GeometryData & geometryData_;

  std::shared_ptr<const orcamodel::Geometry> geom_;

  size_t ctlVecSize_;
  std::shared_ptr<nv::GeometryNV> nvgeom_;
  std::shared_ptr<nv::VariablesNV> nvvars_;
  std::shared_ptr<nv::ErrorCovarianceNV> nverrorcov_;
  std::string covtyp_;

  int nx_;
  int ny_;

//  void print(std::ostream &) const override;

/// inner Geometry Data for next block
//  const oops::GeometryData & innerGeometryData_;
/// inner variables for next block
//  const oops::Variables innerVars_;
};

}  // namespace orcajedifilter
}  // namespace saber
