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
#include "orca-jedi/geometry/GeometryParameterTraitsFieldDType.h"
#include "orca-jedi/utilities/Types.h"

namespace orcamodel {

class NemoFieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NemoFieldParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> name {"name", this};
  oops::RequiredParameter<std::string> nemoName {"nemo field name", this};
  oops::RequiredParameter<std::string> modelSpace {"model space", this};
  oops::Parameter<std::string> variableType {"variable type",
    "type of variable (default is 'background' other options are 'background error variance' and"
    "'background error standard deviation' both are included for clarity, but both variables are"
    " read from the error file",
    "background",
    this};
  oops::Parameter<FieldDType> fieldPrecision{"field precision",
    "Precision to store atlas fields (float (default) or double).",
    FieldDType::Float,
    this};
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
  oops::OptionalParameter<bool> extraFieldsInit{"extrafields initialisation", this};
};

}  //  namespace orcamodel
