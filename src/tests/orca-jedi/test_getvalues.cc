#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/getvalues/GetValues.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE ("test basic getvalues") {
    
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

  eckit::LocalConfiguration getvalues_conf;  
  std::vector<float> lons{0, 120, 270};
  std::vector<float> lats{88, 0, 30};
  std::vector<util::DateTime> times{
    util::DateTime("2018-04-15T00:00:00Z"),
    util::DateTime("2018-04-15T00:00:00Z"),
    util::DateTime("2018-04-15T00:00:00Z")
  };
  std::shared_ptr<const ioda::Distribution> distribution;

  ufo::Locations locations(lons, lats, times, distribution);
  GetValues getvalues(geometry, locations, getvalues_conf);

  // create a state from the test data 
  eckit::LocalConfiguration state_config;  
  std::vector<std::string> state_variables {"sea_ice_area_fraction"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2018-04-15T00:00:00Z");
  state_config.set("nemo field file", "../testinput/orca2_t_nemo.nc");
  state_config.set("variance field file", "../testinput/orca2_t_bkg_var.nc");
  State state(geometry, state_config);

  // create geovals from the locations 
  std::vector<size_t> nlevs = {1};
  ufo::GeoVaLs geovals(locations, state.variables(), std::vector<size_t>{1, 1});

  util::DateTime dt_begin("2018-04-14T00:00:00Z");
  util::DateTime dt_end("2018-04-16T00:00:00Z");
  getvalues.fillGeoVaLs(state, dt_begin, dt_end, geovals);

  // test fillGeoVaLs
  std::vector<double> vals(geovals.nlocs());
  geovals.get(vals, "sea_ice_area_fraction");

  std::vector<double> testvals = {1, 0, 0};

  for (size_t i=0; i<vals.size(); ++i) {
    EXPECT_EQUAL(vals[i], testvals[i]);
  }

}

}  // namespace test
}  // namespace orcamodel

int main( int argc, char** argv ) {
    return orcamodel::test::run( argc, argv );
}
