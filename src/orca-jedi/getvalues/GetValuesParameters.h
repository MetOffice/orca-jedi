/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace orcamodel {

class OrcaAtlasInterpolatorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaAtlasInterpolatorParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> type{"type", this};
  oops::RequiredParameter<std::string> non_linear{"non_linear", this};
};

class OrcaGetValuesParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaGetValuesParameters, oops::Parameters)

 public:
  oops::RequiredParameter<OrcaAtlasInterpolatorParameters> atlasInterpolator{
    "atlas-interpolator", this};
};

}  //  namespace orcamodel
