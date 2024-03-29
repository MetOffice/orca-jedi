/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "atlas/library/Library.h"

#include "orca-jedi/increment/Increment.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test create increment") {
  EXPECT(eckit::system::Library::exists("atlas-orca"));

  eckit::LocalConfiguration config;
  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
    .set("nemo field name", "sic_tot_var")
    .set("model space", "surface");
  nemo_var_mappings[2].set("name", "sea_surface_foundation_temperature")
    .set("nemo field name", "votemper")
    .set("model space", "surface");
  nemo_var_mappings[3].set("name", "sea_water_potential_temperature")
    .set("nemo field name", "votemper")
    .set("model space", "volume");
  config.set("nemo variables", nemo_var_mappings);
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 10);
  Geometry geometry(config, eckit::mpi::comm());

  const std::vector<int> channels{};
  std::vector<std::string> varnames {"sea_ice_area_fraction",
    "sea_water_potential_temperature"};
  oops::Variables oops_vars(varnames, channels);

  util::DateTime datetime;

  EXPECT_THROWS_AS(Increment increment(geometry, oops_vars, datetime),
      eckit::NotImplemented);
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
