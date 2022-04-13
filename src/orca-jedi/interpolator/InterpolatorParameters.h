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
  oops::OptionalParameter<std::string> non_linear{"non_linear", this};
  oops::OptionalParameter<double>
    max_fraction_elems_to_try{"max_fraction_elems_to_try", this};
};

class OrcaInterpolatorParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaInterpolatorParameters, oops::Parameters)

 public:
  oops::RequiredParameter<OrcaAtlasInterpolatorParameters> atlasInterpolator{
    "atlas-interpolator", this};
};

}  //  namespace orcamodel
