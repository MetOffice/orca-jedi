#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE ("test basic state") {
    
  eckit::LocalConfiguration config;  
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 2);

  eckit::LocalConfiguration nemo_var_mapping;  
  nemo_var_mapping.set("sea_ice_area_fraction", "iiceconc");
  nemo_var_mapping.set("sea_surface_foundation_temperature", "votemper");
  nemo_var_mapping.set("sea_water_potential_temperature", "votemper");
  config.set("nemo names", nemo_var_mapping);
  config.set("variance names", std::vector<std::string>{"sic_tot_var"});
  Geometry geometry(config, eckit::mpi::comm());
  const std::vector<int> channels{};

  eckit::LocalConfiguration state_config;  
  std::vector<std::string> state_variables {"sea_ice_area_fraction"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2018-04-15T00:00:00Z");

  SECTION ("test constructor") {
    oops::Variables oops_vars(state_variables, channels);
    util::DateTime datetime("2018-04-15T00:00:00Z");
    State state(geometry, oops_vars, datetime);
  }


  SECTION ("test constructor from config") {
    state_config.set("nemo field file", "../testinput/orca2_t_nemo.nc");
    state_config.set("variance field file", "../testinput/orca2_t_bkg_var.nc");
    State state(geometry, state_config);
    bool has_missing = state.stateFields()["iiceconc"].metadata().has("missing_value");
    EXPECT_EQUAL(true, has_missing);
  }

  SECTION ("test constructor from config analytic_init") {
    state_config.set("nemo field file", "../testinput/orca2_t_nemo.nc");
    state_config.set("analytic_init", "zeroed state");
    State state(geometry, state_config);
  }

  EXPECT(true);

}

}  // namespace test
}  // namespace orcamodel

int main( int argc, char** argv ) {
    return orcamodel::test::run( argc, argv );
}
