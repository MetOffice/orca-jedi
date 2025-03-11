/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include "oops/base/Variables.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/OptionalParameter.h"


class VariableChangeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters, oops::Parameters)
 public:
  oops::OptionalParameter<oops::Variables> inputVariables{"input variables", this};
  oops::OptionalParameter<oops::Variables> outputVariables{"output variables", this};
  // No linear variable change. No additional parameters
};

