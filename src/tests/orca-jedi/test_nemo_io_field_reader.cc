/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test non-existent file throws error") {
  eckit::PathName test_data_path("NOTAFILE.not_nc");
  EXPECT_THROWS_AS(NemoFieldReader field_reader(test_data_path),
      eckit::BadValue);
}

CASE("test corrupted file throws error") {
  eckit::PathName test_data_path("../Data/hofx3d_nc_sst.yaml");
  EXPECT_THROWS_AS(NemoFieldReader field_reader(test_data_path),
      eckit::BadValue);
}

CASE("test opening bad test file ") {
  eckit::PathName test_data_path("../Data/orca2_t_coords.nc");
  EXPECT_THROWS_AS(NemoFieldReader field_reader(test_data_path),
      eckit::BadValue);
}

CASE("test opening the test file ") {
  eckit::PathName test_data_path("../Data/simple_nemo.nc");
  NemoFieldReader field_reader(test_data_path);
  SECTION("reading missing dimension throws error") {
    EXPECT_THROWS_AS(field_reader.read_dim_size("NOTADIMENSION"),
        eckit::BadValue);
  }
}

CASE("test reading the latitudes and longitudes arrays ") {
  eckit::PathName test_data_path("../Data/simple_nemo.nc");

  NemoFieldReader field_reader(test_data_path);
  std::vector<atlas::PointXY> ob_locs = field_reader.read_locs();

  EXPECT_EQUAL(ob_locs.size(), 6);

  EXPECT_EQUAL(ob_locs[0].x(), 10);
  EXPECT_EQUAL(ob_locs[0].y(), 0);
  EXPECT_EQUAL(ob_locs[5].x(), 12);
  EXPECT_EQUAL(ob_locs[5].y(), 1);
}

CASE("test get_nearest_datetime_index returns correct index") {
  eckit::PathName test_data_path("../Data/simple_nemo.nc");

  NemoFieldReader field_reader(test_data_path);
  util::DateTime test_datetime{"1970-01-01T00:00:00Z"};
  size_t index = field_reader.get_nearest_datetime_index(test_datetime);

  EXPECT_EQUAL(index, 0);
}

CASE("test reading field _FillValue") {
  eckit::PathName test_data_path("../Data/simple_nemo.nc");

  NemoFieldReader field_reader(test_data_path);

  auto missing_value = field_reader.read_fillvalue<double>("iiceconc");
  EXPECT_EQUAL(missing_value, -32768.);

  auto default_missing_value = field_reader.read_fillvalue<double>("nav_lat");
  EXPECT_EQUAL(default_missing_value, std::numeric_limits<double>::lowest());
}

CASE("test read_var_slice reads vector") {
  eckit::PathName test_data_path("../Data/simple_nemo.nc");

  NemoFieldReader field_reader(test_data_path);
  std::vector<double> data = field_reader.read_var_slice<double>("iiceconc", 1, 0);

  EXPECT_EQUAL(data[0], 121);
  EXPECT_EQUAL(data[5], 171);
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
