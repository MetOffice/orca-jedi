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
#include "oops/util/parameters/OptionalParameter.h"

namespace orcamodel {

class NemoFieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NemoFieldParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> name {"name", this};
  oops::RequiredParameter<std::string> nemoName {"nemo field name", this};
  oops::RequiredParameter<std::string> modelSpace {"model space", this};
  oops::OptionalParameter<std::string> variableType {"variable type", this};
};

class OrcaGeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(OrcaGeometryParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::vector<NemoFieldParameters>> nemoFields
    {"nemo variables", this};
  oops::RequiredParameter<std::string> gridName
    {"grid name", this};
  oops::RequiredParameter<int> nLevels {"number levels", this};
  oops::Parameter<int> sourceMeshHalo {"source mesh halo",
    "Size of the MPI halo when using a domain-distributed geometry."
      " The default is 0 (no MPI halo)",
    0,
    this};
  oops::Parameter<std::string> partitioner {
    "partitioner",
      "Name of the atlas partitioner to use to MPI distribute the model data"
        " The default will not distribute the data ('serial').",
      "serial",
      this};
};

}  //  namespace orcamodel
