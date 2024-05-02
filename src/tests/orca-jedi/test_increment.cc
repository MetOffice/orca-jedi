/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <iostream>

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

CASE("test basic increment") {
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
  std::vector<std::string> varnames2 {"sea_ice_area_fraction",
    "sea_water_potential_temperature"};
  oops::Variables oops_vars2(varnames2, channels);

  std::vector<std::string> varnames {"sea_ice_area_fraction"};
  oops::Variables oops_vars(varnames, channels);

  util::DateTime datetime("2021-06-30T00:00:00Z");

  SECTION("test constructor") {
    Increment increment(geometry, oops_vars2, datetime);
  }

  SECTION("test setting increment value") {
    std::cout << "DJL test setting increment value a" << std::endl;
    Increment increment(geometry, oops_vars, datetime);
    std::cout << "DJL test setting increment value b" << std::endl;
    increment.ones();
    increment.print(std::cout);
    std::cout << "DJL test setting increment value c" << std::endl;
    EXPECT_EQUAL(increment.norm(), 1);
    std::cout << "DJL test setting increment value d" << std::endl;
    increment.zero();
    increment.print(std::cout);
    EXPECT_EQUAL(increment.norm(), 0);
    std::cout << "DJL test setting increment value e" << std::endl;
  }

  SECTION("test dirac") {
    std::cout << "DJL test dirac" << std::endl;
    eckit::LocalConfiguration dirac_config;
    std::vector<int> ix = {20, 30};
    std::vector<int> iy = {10, 40};
    std::vector<int> iz = {1, 3};
    dirac_config.set("ixdir", ix);
    dirac_config.set("iydir", iy);
    dirac_config.set("izdir", iz);

    Increment increment(geometry, oops_vars, datetime);
    increment.dirac(dirac_config);
    increment.print(std::cout);
    
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
