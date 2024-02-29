/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/state/StateIOUtils.h"

#include <sstream>
#include <vector>
#include <map>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "atlas/field.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/nemo_io/NemoFieldWriter.h"

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
    ReadServer nemo_reader(geom.timer(), nemo_field_path, geom.mesh());

    // Read fields from Nemo field file
    // field names in the atlas fieldset are assumed to match their names in
    // the field file
    const size_t time_indx = nemo_reader.get_nearest_datetime_index(valid_date);
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
        const auto populate = [&](auto typeVal) {
          using T = decltype(typeVal);
          populateField<T>(nemoName, varCoordTypeMap[fieldName],
                                time_indx, nemo_reader, field);
        };
        ApplyForFieldType(populate,
                          geom.fieldPrecision(fieldName),
                          eckit::BadParameter("State(ORCA)::readFieldsFromFile '"
                            + nemoName + "' field type not recognised"));
        // Add a halo exchange following read to fill out halo points
        geom.functionSpace().haloExchange(field);
        geom.log_status();
      }
    }

    oops::Log::trace() << "orcamodel::readFieldsFromFile:: readFieldsFromFile "
                       << "done" << std::endl;
}

/// \brief Populate a single atlas field using the read server.
/// \param nemo_name The netCDF name of the variable to read.
/// \param coord_type The type of coordinate (e.g "vertical" for 1D data).
/// \param time_indx The time index in the file.
/// \param nemo_reader The read server managing IO with the file.
template<class T> void populateField(
  const std::string & nemo_name,
  const std::string & coord_type,
  size_t time_indx,
  ReadServer & nemo_reader,
  atlas::Field & field) {
    atlas::array::ArrayView<T, 2> field_view =
        atlas::array::make_view<T, 2>(field);
    if (coord_type == "vertical") {
      nemo_reader.read_vertical_var<T>(nemo_name, field_view);
    } else {
      nemo_reader.read_var<T>(nemo_name, time_indx, field_view);
    }
    T missing_value = nemo_reader.read_fillvalue<T>(nemo_name);
    field.metadata().set("missing_value", missing_value);
    field.metadata().set("missing_value_type", "approximately-equals");
    field.metadata().set("missing_value_epsilon", NEMO_FILL_TOL);
}
template void populateField<double>(
  const std::string & nemo_name,
  const std::string & coord_type,
  size_t time_indx,
  ReadServer & nemo_reader,
  atlas::Field & field);
template void populateField<float>(
  const std::string & nemo_name,
  const std::string & coord_type,
  size_t time_indx,
  ReadServer & nemo_reader,
  atlas::Field & field);

void writeFieldsToFile(
  const OrcaStateParameters & params,
  const Geometry & geom,
  const util::DateTime & valid_date,
  const atlas::FieldSet & fs) {
    oops::Log::trace() << "orcamodel::writeFieldsToFile:: start for valid_date"
                       << " " << valid_date << std::endl;

    std::string output_filename =
      params.outputNemoFieldFile.value().value_or("");
    if (output_filename == "")
      throw eckit::BadValue(std::string("orcamodel::writeFieldsToFile:: ")
          + "file name not specified", Here());

    std::map<std::string, std::string> varCoordTypeMap;
    {
      const oops::Variables vars = geom.variables();
      const std::vector<std::string> coordSpaces =
        geom.variableNemoSpaces(vars);
      for (int i=0; i < vars.size(); ++i)
        varCoordTypeMap[vars[i]] = coordSpaces[i];
    }

    auto nemo_field_path = eckit::PathName(output_filename);
    oops::Log::debug() << "orcamodel::writeFieldsToFile:: nemo_field_path "
                       << nemo_field_path << std::endl;
    std::vector<util::DateTime> datetimes = {valid_date};
    std::vector<double> levels((*fs.begin()).shape(1), 0);
    for (size_t iLev = 0; iLev < levels.size(); ++iLev) { levels[iLev] = iLev; }

    auto writeRankFields = [&](){
      NemoFieldWriter field_writer(nemo_field_path, geom.mesh(), datetimes,
          levels);
      for (atlas::Field field : fs) {
        std::string fieldName = field.name();
        std::string nemoName = geom.nemo_var_name(fieldName);
        auto field_view = atlas::array::make_view<double, 2>(field);
        if (varCoordTypeMap[fieldName] == "surface") {
          field_writer.write_surf_var(nemoName, field_view, 0);
        } else {
          field_writer.write_vol_var(nemoName, field_view, 0);
        }
      }
    };

    // Write from rank 0 unless the data is distributed. If the date is
    // distributed sequentially write the data to the output file.
    if (geom.distributionType() == "serial") {
        if (geom.getComm().rank() == 0) writeRankFields();
    } else {
      geom.getComm().barrier();
      size_t rank = 0;
      while (rank < geom.getComm().size()) {
        if (rank == geom.getComm().rank()) {
          writeRankFields();
        }
        rank++;
      }
    }
    geom.getComm().barrier();
}

}  // namespace orcamodel
