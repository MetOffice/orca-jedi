/*
 * (C) British Crown Copyright 2020-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include<sstream>

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateParameters.h"
#include "orca-jedi/interpolator/Interpolator.h"
#include "orca-jedi/interpolator/InterpolatorParameters.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic interpolator") {
  eckit::LocalConfiguration config;
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 2);

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
  Geometry geometry(config, eckit::mpi::comm());

  eckit::LocalConfiguration interp_conf;
  interp_conf.set("type", "finite-element");
  interp_conf.set("non_linear", "missing-if-all-missing");
  eckit::LocalConfiguration interpolator_conf;
  interpolator_conf.set("atlas-interpolator", interp_conf);

  OrcaInterpolatorParameters params;
  params.validateAndDeserialize(interpolator_conf);

  // lons{0, 120, 270};
  // lats{88, 0, 30};
  std::vector<double> locations({88, 0, 0, 120, 30, 270});

  // create a state from the test data
  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {"sea_ice_area_fraction", "sea_surface_foundation_temperature"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2018-04-15T00:00:00Z");
  state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
  state_config.set("variance field file", "../Data/orca2_t_bkg_var.nc");
  OrcaStateParameters stateParams;
  stateParams.validateAndDeserialize(state_config);
  State state(geometry, stateParams);

  SECTION("test get values fails with no locations") {
    EXPECT_THROWS_AS(
      Interpolator interpolator(interpolator_conf, geometry, {}),
      eckit::BadValue);
  }

  Interpolator interpolator(interpolator_conf, geometry, locations);

  SECTION("test interpolator.apply fails missing variable") {
    oops::Variables variables({"NOTAVARIABLE"});
    std::vector<double> vals(3);
    EXPECT_THROWS_AS(interpolator.apply(variables, state, vals),
        eckit::BadParameter);
  }

  SECTION("test interpolator.apply") {

    std::vector<double> vals(locations.size() / 2);
    interpolator.apply(oops::Variables({"sea_ice_area_fraction", "sea_surface_foundation_temperature"}), state, vals);

    double missing_value = util::missingValue(vals[0]);
    std::vector<double> testvals = {1, missing_value, 0, 18, 18, 18};

    EXPECT_EQUAL(vals[0], testvals[0]);
    EXPECT_EQUAL(vals[1], testvals[1]);
    EXPECT_EQUAL(vals[2], testvals[2]);
    //EXPECT_EQUAL(vals[3], testvals[3]);
    //EXPECT_EQUAL(vals[4], testvals[4]);
    //EXPECT_EQUAL(vals[5], testvals[5]);
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run( argc, argv );
}
