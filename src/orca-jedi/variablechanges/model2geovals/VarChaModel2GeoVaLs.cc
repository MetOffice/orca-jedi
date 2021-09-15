/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/variablechanges/model2geovals/VarChaModel2GeoVaLs.h"

namespace orcamodel {
// -------------------------------------------------------------------------------------------------
static oops::VariableChangeMaker<OrcaModelTraits, VarChaModel2GeoVaLs>
       makerVarChaModel2GeoVaLs_("Model2GeoVaLs");
static oops::VariableChangeMaker<OrcaModelTraits, VarChaModel2GeoVaLs> makerVarChaDefault_("default");
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::VarChaModel2GeoVaLs(const Geometry & geom, const eckit::Configuration & conf) :
  geom_(new Geometry(geom))
{
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaModel2GeoVaLs::~VarChaModel2GeoVaLs() {
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVar(const State & xin, State & xout) const {
  oops::Log::trace() << classname() << " changeVar start" << std::endl;
  xout = xin;
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::changeVarInverse(const State & xin, State & xout) const {
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  xout = xin;
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaModel2GeoVaLs::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace orcamodel
