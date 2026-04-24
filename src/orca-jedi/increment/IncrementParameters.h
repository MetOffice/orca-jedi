/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace orcamodel {

class OrcaIncrementParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaIncrementParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> nemoFieldFile{"output path", this};
  oops::RequiredParameter<util::DateTime> date{"date", this};
};

class OrcaDiracParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaDiracParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::vector<int>> ixdir{"x indices", this};
  oops::RequiredParameter<std::vector<int>> iydir{"y indices", this};
  oops::RequiredParameter<std::vector<int>> izdir{"z indices", this};
};
}  //  namespace orcamodel
