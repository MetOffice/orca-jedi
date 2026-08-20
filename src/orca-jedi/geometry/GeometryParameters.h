/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"  // IWYU pragma: keep
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "orca-jedi/geometry/GeometryParameterTraitsFieldDType.h"  // IWYU pragma: keep
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

/// \brief Optional parameters for reading a land-sea mask ancillary.
///
/// When specified in the geometry configuration, a field is read from the
/// given NetCDF file and used to populate the 3D volume_mask extra field. Two
/// modes are supported via "type":
///   - "missing value" (default): the field is a data field and its missing
///     values define the land points.
///   - "mask_field": the field is an explicit land-sea mask; entries approximately
///     equal to "land value" (default 0) are land, all others are ocean.
/// This makes the mask available to any code that uses the geometry
/// (e.g. State resolution-change constructor).
class LandSeaMaskParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LandSeaMaskParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> filepath {"filepath", this};
  oops::RequiredParameter<std::string> variable {"variable", this};
  oops::Parameter<std::string> type {"type",
    "How to derive the mask from the variable: 'missing value' (mask where the"
    " field is missing, the default) or 'mask_field' (use the field values"
    " directly).",
    "missing value", this};
  oops::Parameter<double> landValue {"land value",
    "For type 'mask_field', field values approximately equal to this are treated"
    " as land (masked). Default 0.",
    0.0, this};
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
  oops::OptionalParameter<bool> extraFieldsInit{
      "initialise extra fields", this};
  oops::OptionalParameter<LandSeaMaskParameters> landSeaMask{
      "land sea mask", this};

  oops::Parameter<bool> parallelOutput{"parallel output",
    "Write NEMO output files using the parallel (MPI + parallel-netCDF) I/O path"
      " instead of gathering every field onto the root rank. Requires an MPI run."
      " The default is false (gather-on-root).",
    false,
    this};
  oops::Parameter<int> outputIoRanks{"output io ranks",
    "Number of ranks to use for parallel output when 'parallel output' is true."
      " A value of 0 (the default) uses every rank in the communicator. This"
      " should typically be tuned towards the number of filesystem stripes.",
    0,
    this};
  oops::Parameter<bool> parallelInput{"parallel input",
    "Read NEMO input files using the parallel (MPI + parallel-netCDF) I/O path"
      " instead of reading on the root rank and broadcasting. Requires an MPI"
      " run and a netCDF-4 (HDF5) input file. The default is false"
      " (read-on-root).",
    false,
    this};
  oops::Parameter<int> inputIoRanks{"input io ranks",
    "Number of ranks to use for parallel input when 'parallel input' is true."
      " A value of 0 (the default) uses every rank in the communicator.",
    0,
    this};
};

}  //  namespace orcamodel
