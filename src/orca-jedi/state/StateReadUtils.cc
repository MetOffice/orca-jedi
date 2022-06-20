/*
 * (C) British Crown Copyright 2020-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "orca-jedi/state/StateReadUtils.h"

#include <sstream>
#include <vector>
#include <map>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"

#define NEMO_FILL_TOL 1e-6

namespace orcamodel {

void readFieldsFromFile(
  const OrcaStateParameters & params,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const std::string & variable_type,
  atlas::FieldSet & fs) {
    oops::Log::trace() << "orcamodel::readFieldsFromFile:: start for valid_date"
                       << " " << valid_date << std::endl;

    // Open Nemo Feedback file

    oops::Log::debug() << "orcamodel::readFieldsFromFile:: parameters "
                       << params
                       << std::endl;

    std::string nemo_file_name;
    if (variable_type == "background") {
      nemo_file_name = params.nemoFieldFile.value();
    }
    if (variable_type == "background variance") {
      nemo_file_name = params.varianceFieldFile.value().value_or("");
    }

    auto nemo_field_path = eckit::PathName(nemo_file_name);
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: nemo_field_path "
                       << nemo_field_path << std::endl;
    NemoFieldReader nemo_file(nemo_field_path);

    // Read fields from Nemo field file
    // field names in the atlas fieldset are assumed to match their names in
    // the field file
    const size_t time_indx = nemo_file.get_nearest_datetime_index(valid_date);
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: time_indx "
                       << time_indx << std::endl;

    std::map<std::string, std::string> varCoordTypeMap;
    {
      const oops::Variables vars = geom.variables();
      const std::vector<std::string> coordSpaces =
        geom.variableNemoSpaces(vars);
      for (int i=0; i < vars.size(); ++i)
        varCoordTypeMap[vars[i]] = coordSpaces[i];
    }
    for (atlas::Field field : fs) {
      std::string fieldName = field.name();
      std::string nemoName = geom.nemo_var_name(fieldName);
      oops::Log::debug() << "orcamodel::readFieldsFromFile:: "
                         << "geom.variable_in_variable_type(\""
                         << fieldName << "\", \"" << variable_type << "\") "
                         << geom.variable_in_variable_type(fieldName,
                              variable_type)
                         << std::endl;
      if (geom.variable_in_variable_type(fieldName, variable_type)) {
        auto field_view = atlas::array::make_view<double, 2>(field);
        if (varCoordTypeMap[fieldName] == "surface") {
          nemo_file.read_surf_var(nemoName, time_indx, field_view);
        } else if (varCoordTypeMap[fieldName] == "vertical") {
          nemo_file.read_vertical_var(nemoName, field_view);
        } else {
          nemo_file.read_volume_var(nemoName, time_indx, field_view);
        }
        auto missing_value = nemo_file.read_fillvalue<double>(nemoName);
        field.metadata().set("missing_value", missing_value);
        field.metadata().set("missing_value_type", "approximately-equals");
        field.metadata().set("missing_value_epsilon", NEMO_FILL_TOL);
      }
    }

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: readFieldsFromFile "
                       << "done" << std::endl;
}

}  // namespace orcamodel
