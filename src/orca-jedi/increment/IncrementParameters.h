/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/base/Variables.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
// include "oops/util/parameters/OptionalParameter.h"

namespace orcamodel {

class OrcaIncrementParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaIncrementParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> nemoFieldFile{"filepath", this};  // DJL reconsider names
  oops::RequiredParameter<util::DateTime> date{"date", this};   // DJL needed?
//  oops::RequiredParameter<oops::Variables> stateVariables{"state variables",
//    this};
};

// DJL Add something for the dirac test parameters  ixdir, iydir, izdir ??

}  //  namespace orcamodel
