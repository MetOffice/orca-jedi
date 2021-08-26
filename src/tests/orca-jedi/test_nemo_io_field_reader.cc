#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "eckit/testing/Test.h"

#include "orca-jedi/nemo_io/NemoFieldReader.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE ("test opening the test file ") {

  eckit::PathName test_data_path("../testinput/simple_nemo.nc");

  NemoFieldReader field_reader( test_data_path );

}

CASE ("test reading the latitudes and longitudes arrays ") {

  eckit::PathName test_data_path("../testinput/simple_nemo.nc");

  NemoFieldReader field_reader( test_data_path );
  std::vector<atlas::PointXY> ob_locs = field_reader.read_locs();

  EXPECT_EQUAL(ob_locs.size(), 6);

  EXPECT_EQUAL(ob_locs[0].x(), 10);
  EXPECT_EQUAL(ob_locs[0].y(), 0);
  EXPECT_EQUAL(ob_locs[5].x(), 12);
  EXPECT_EQUAL(ob_locs[5].y(), 1);

}

CASE ("test get_nearest_datetime_index returns correct index") {

  eckit::PathName test_data_path("../testinput/simple_nemo.nc");

  NemoFieldReader field_reader( test_data_path );
  util::DateTime test_datetime{"1970-01-01T00:00:00Z"};
  size_t index = field_reader.get_nearest_datetime_index( test_datetime );

  EXPECT_EQUAL( index, 0 );

}

CASE ("test read_surf_var reads vector") {

  eckit::PathName test_data_path("../testinput/simple_nemo.nc");

  NemoFieldReader field_reader( test_data_path );
  std::vector<double> data = field_reader.read_surf_var( "iiceconc" );

  EXPECT_EQUAL( data[0], 120 );
  EXPECT_EQUAL( data[5], 170 );

}

CASE ("test read_surf_var reads field array view") {

  eckit::PathName test_data_path("../testinput/simple_nemo.nc");

  atlas::idx_t num_points = 6;
  auto array_shape = atlas::array::ArrayShape{6};
  auto array_datatype = atlas::array::DataType::real64();
  std::unique_ptr<atlas::array::Array> test_array(
    atlas::array::Array::create(array_datatype, array_shape) );

  atlas::array::ArrayView<double, 1> arrayView =
    atlas::array::make_view<double, 1>( *test_array );

  NemoFieldReader field_reader( test_data_path );
  field_reader.read_surf_var( "iiceconc", arrayView );

  EXPECT_EQUAL( arrayView(0), 120 );
  EXPECT_EQUAL( arrayView(5), 170 );

}


CASE ("example test") {
    struct Section {
        std::string gridname;
        size_t size;
    };

    std::vector<Section> sections{
        {"ORCA2_T", 27118}, {"eORCA1_T", 120184},
    };
    for ( auto& section : sections ) {
        std::string gridname = section.gridname;
        SECTION( gridname ) {

            size_t n = section.size;
            EXPECT_EQUAL( n, section.size );

        }
    }
}
}  // namespace test
}  // namespace orcamodel

int main( int argc, char** argv ) {
    return eckit::testing::run_tests( argc, argv );
}
