/*
 * (C) British Crown Copyright 2023 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/base/Variables.h"
#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/OptionalParameter.h"


class VariableChangeParameters :
    public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParameters,
                           oops::VariableChangeParametersBase)
 public:
  // No linear variable change. No additional parameters
};

