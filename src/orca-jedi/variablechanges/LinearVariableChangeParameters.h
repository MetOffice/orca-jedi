/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#pragma once

#include "oops/base/Variables.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/OptionalParameter.h"


class LinearVariableChangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters, oops::Parameters)
 public:
  oops::OptionalParameter<oops::Variables> inputVariables{"input variables", this};
  oops::OptionalParameter<oops::Variables> outputVariables{"output variables", this};
};

