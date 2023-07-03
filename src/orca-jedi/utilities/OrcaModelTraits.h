/*
 * (C) British Crown Copyright 2017-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <string>

#include "oops/generic/AtlasInterpolator.h"

#include "orca-jedi/errorcovariance/ErrorCovariance.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"

#include "orca-jedi/model/ModelBias.h"
#include "orca-jedi/model/ModelBiasIncrement.h"
#include "orca-jedi/model/ModelBiasCovariance.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/variablechanges/VariableChange.h"

namespace orcamodel {

struct OrcaModelTraits {
  static std::string name() {return "OrcaModel";}
  static std::string nameCovar() {return "ORCAstatic";}
  static std::string nameCovar4D() {return "ORCAstatic";}

  typedef orcamodel::ErrorCovariance           Covariance;
  typedef orcamodel::Geometry                  Geometry;

  typedef oops::AtlasInterpolator              LocalInterpolator;
  typedef orcamodel::Increment                 Increment;

  typedef orcamodel::ModelBias                 ModelAuxControl;
  typedef orcamodel::ModelBiasIncrement        ModelAuxIncrement;
  typedef orcamodel::ModelBiasCovariance       ModelAuxCovariance;
  typedef orcamodel::State                     State;
  typedef orcamodel::VariableChange            VariableChange;
};

}  // namespace orcamodel
