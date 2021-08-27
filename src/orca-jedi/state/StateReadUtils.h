
#pragma once

#include "eckit/config/Configuration.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

namespace orcamodel {

void readFieldsFromFile(
  const eckit::Configuration & conf,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const std::string & variable_type,
  atlas::FieldSet & fs);

}
