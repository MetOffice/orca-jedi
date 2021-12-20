/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/


#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/base/Variables.h"
#include "oops/base/LinearVariableChangeParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/OptionalParameter.h"


class LinearVariableChangeParameters :
    public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(LinearVariableChangeParameters,
                           oops::LinearVariableChangeParametersBase)
 public:
  // No linear variable change. No additional parameters
};

