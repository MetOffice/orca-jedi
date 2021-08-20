#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE ("test basic geometry") {
    
  EXPECT( eckit::system::Library::exists( "atlas-orca" ) );

  eckit::LocalConfiguration config;  
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 10);
  eckit::LocalConfiguration nemo_var_mapping;  
  nemo_var_mapping.set("sea_ice_area_fraction", "iiceconc");
  nemo_var_mapping.set("sea_surface_foundation_temperature", "votemper");
  nemo_var_mapping.set("sea_water_potential_temperature", "votemper");
  config.set("nemo names", nemo_var_mapping);
  Geometry geometry(config, eckit::mpi::comm());

  SECTION( "test geometry variable sizes" ) {

    const std::vector<int> channels{};
    std::vector<std::string> varnames {"sea_ice_area_fraction", "sea_water_potential_temperature"};
    oops::Variables oops_vars(varnames, channels);
    auto varsizes = geometry.variableSizes( oops_vars );
    EXPECT_EQUAL(varsizes.size(), 2);
    EXPECT_EQUAL(varsizes[0], 1);
    EXPECT_EQUAL(varsizes[1], 10);

  }
}

}  // namespace test
}  // namespace orcamodel

int main( int argc, char** argv ) {
    return orcamodel::test::run( argc, argv );
}
