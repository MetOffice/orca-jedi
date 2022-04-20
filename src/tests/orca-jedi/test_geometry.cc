/*
 * (C) British Crown Copyright 2020-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic geometry") {
  EXPECT(eckit::system::Library::exists("atlas-orca"));

  eckit::LocalConfiguration config;
  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
    .set("nemo field name", "sic_tot_var")
    .set("model space", "surface")
    .set("variable type", "background variance");
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

  SECTION("test geometry variable names") {
    EXPECT_THROWS_AS(geometry.nemo_var_name("NOTAVARIABLE"), eckit::BadValue);
    EXPECT(geometry.variable_in_variable_type("sea_ice_area_fraction",
                                              "background"));
    EXPECT(geometry.variable_in_variable_type("sea_ice_area_fraction_error",
                                              "background variance"));
    EXPECT(!geometry.variable_in_variable_type("sea_ice_area_fraction_error",
                                               "background"));
    EXPECT(!geometry.variable_in_variable_type("sea_ice_area_fraction",
                                               "background variance"));
  }

  SECTION("test geometry variable sizes") {
    const std::vector<int> channels{};
    std::vector<std::string> varnames {"sea_ice_area_fraction",
      "sea_water_potential_temperature"};
    oops::Variables oops_vars(varnames, channels);
    auto varsizes = geometry.variableSizes(oops_vars);
    EXPECT_EQUAL(varsizes.size(), 2);
    EXPECT_EQUAL(varsizes[0], 1);
    EXPECT_EQUAL(varsizes[1], 10);
    oops::Variables not_vars({"NOTAVARIBLE"}, channels);
    EXPECT_THROWS_AS(geometry.variableSizes(not_vars), eckit::BadValue);
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
