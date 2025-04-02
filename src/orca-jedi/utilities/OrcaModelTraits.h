/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>

#include "orca-jedi/errorcovariance/ErrorCovariance.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/interpolator/Interpolator.h"
#include "orca-jedi/increment/Increment.h"

#include "orca-jedi/model/ModelBias.h"
#include "orca-jedi/model/ModelBiasIncrement.h"
#include "orca-jedi/model/ModelBiasCovariance.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/utilities/ModelData.h"
#include "orca-jedi/variablechanges/VariableChange.h"
#include "orca-jedi/variablechanges/LinearVariableChange.h"

namespace orcamodel {

struct OrcaModelTraits {
  static std::string name() {return "OrcaModel";}
  static std::string nameCovar() {return "ORCAstatic";}
  static std::string nameCovar4D() {return "ORCAstatic";}

  typedef orcamodel::ErrorCovariance           Covariance;
  typedef orcamodel::Geometry                  Geometry;

  typedef orcamodel::Interpolator              LocalInterpolator;
  typedef orcamodel::Increment                 Increment;

  typedef orcamodel::ModelBias                 ModelAuxControl;
  typedef orcamodel::ModelBiasIncrement        ModelAuxIncrement;
  typedef orcamodel::ModelBiasCovariance       ModelAuxCovariance;
  typedef orcamodel::State                     State;
  typedef orcamodel::VariableChange            VariableChange;
  typedef orcamodel::LinearVariableChange      LinearVariableChange;
  typedef orcamodel::ModelData                 ModelData;
};

}  // namespace orcamodel
