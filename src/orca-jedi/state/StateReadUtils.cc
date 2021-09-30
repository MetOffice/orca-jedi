
#include <sstream>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"

namespace orcamodel {

void readFieldsFromFile(
  const eckit::Configuration & conf,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const std::string & variable_type,
  atlas::FieldSet & fs) {

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: start for valid_date "
                       << valid_date << std::endl;

    // Open Nemo Feedback file

    oops::Log::debug() << "orcamodel::readFieldsFromFile:: configuration "
                       << conf
                       << std::endl;

    std::string nemo_file_name;
    if (variable_type == "background" ) {
      nemo_file_name = conf.getString("nemo field file");
    }
    if (variable_type == "background variance" ) {
      nemo_file_name = conf.getString("variance field file");
    }

    auto nemo_field_path = eckit::PathName(nemo_file_name);
    oops::Log::debug() << "orcamodel::readFieldsFromFile:: nemo_field_path "
                       << nemo_field_path << std::endl;
    NemoFieldReader nemo_file(nemo_field_path);

    // Read fields from Nemo field file
    // field names in the atlas fieldset are assumed to match their names in
    // the field file
    size_t time_indx = nemo_file.get_nearest_datetime_index(valid_date);

    for (atlas::Field field : fs) {
      std::string fieldName = field.name();
      oops::Log::debug() << "orcamodel::readFieldsFromFile:: field name = " << fieldName
                         << std::endl;

      oops::Log::debug() << "orcamodel::readFieldsFromFile:: "
                         << "geom.variable_in_variable_type(\"" << fieldName << "\", \"" << variable_type << "\") "
                         << geom.variable_in_variable_type(fieldName, variable_type) << std::endl;
      if (geom.variable_in_variable_type(fieldName, variable_type)) {

        auto field_view = atlas::array::make_view<double, 1>( field );
        nemo_file.read_surf_var(fieldName, time_indx, field_view);

        auto missing_value = nemo_file.read_fillvalue<double>(fieldName);
        field.metadata().set("missing_value", missing_value);
      }
    }

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: readFieldsFromFile done "
                       << std::endl;
};

}
