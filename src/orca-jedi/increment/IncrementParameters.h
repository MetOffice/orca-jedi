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

namespace orcamodel {

class OrcaIncrementParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaIncrementParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> nemoFieldFile{"filepath", this};
  oops::RequiredParameter<util::DateTime> date{"date", this};
};

class OrcaDiracParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaDiracParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::vector<int>> ixdir{"ixdir", this};
  oops::RequiredParameter<std::vector<int>> iydir{"iydir", this};
  oops::RequiredParameter<std::vector<int>> izdir{"izdir", this};
};
}  //  namespace orcamodel
