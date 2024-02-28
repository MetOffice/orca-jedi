/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateParameters.h"
#include "orca-jedi/utilities/Types.h"
#include "orca-jedi/nemo_io/ReadServer.h"

namespace orcamodel {

void readFieldsFromFile(
  const OrcaStateParameters & params,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const std::string & variable_type,
  atlas::FieldSet & fs);
void writeFieldsToFile(
  const OrcaStateParameters & params,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const atlas::FieldSet & fs);
template<class T> void populateField(
  const std::string & nemo_name,
  const std::string & coord_type,
  size_t time_indx,
  ReadServer & nemo_reader,
  atlas::Field & field);

}  // namespace orcamodel
