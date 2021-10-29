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

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateParameters.h"
#include "orca-jedi/getvalues/GetValues.h"
#include "orca-jedi/getvalues/GetValuesParameters.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic getvalues") {
  eckit::LocalConfiguration config;
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 2);

  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
    .set("nemo field name", "iiceconc")
    .set("type", "surface");
  nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
    .set("nemo field name", "sic_tot_var")
    .set("type", "surface");
  nemo_var_mappings[2].set("name", "sea_surface_foundation_temperature")
    .set("nemo field name", "votemper")
    .set("type", "surface");
  nemo_var_mappings[3].set("name", "sea_water_potential_temperature")
    .set("nemo field name", "votemper")
    .set("type", "volume");
  config.set("nemo variables", nemo_var_mappings);
  config.set("variance names", std::vector<std::string>{"sic_tot_var"});
  Geometry geometry(config, eckit::mpi::comm());

  eckit::LocalConfiguration interp_conf;
  interp_conf.set("type", "finite-element");
  interp_conf.set("non_linear", "missing-if-all-missing");
  eckit::LocalConfiguration getvalues_conf;
  getvalues_conf.set("atlas-interpolator", interp_conf);

  OrcaGetValuesParameters params;
  params.validateAndDeserialize(getvalues_conf);

  std::vector<float> lons{0, 120, 270};
  std::vector<float> lats{88, 0, 30};
  std::vector<util::DateTime> times{
    util::DateTime("2018-04-15T00:00:00Z"),
    util::DateTime("2018-04-15T00:00:00Z"),
    util::DateTime("2018-04-15T00:00:00Z")
  };
  std::shared_ptr<const ioda::Distribution> distribution;

  ufo::Locations locations(lons, lats, times, distribution);

  // create a state from the test data
  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {"sea_ice_area_fraction"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2018-04-15T00:00:00Z");
  state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
  state_config.set("variance field file", "../Data/orca2_t_bkg_var.nc");
  OrcaStateParameters stateParams;
  stateParams.validateAndDeserialize(state_config);
  State state(geometry, stateParams);

  util::DateTime dt_begin("2018-04-14T00:00:00Z");
  util::DateTime dt_end("2018-04-16T00:00:00Z");

  SECTION("test get values fails with no locations") {
    ufo::Locations emptyLocations({}, {}, {}, distribution);
    EXPECT_THROWS_AS(GetValues getvalues(geometry, emptyLocations,
        getvalues_conf), eckit::BadValue);
  }

  GetValues getvalues(geometry, locations, getvalues_conf);

  SECTION("test fillGeoVaLs fails missing variable") {
    oops::Variables variables({"NOTAVARIABLE"});
    ufo::GeoVaLs geovals(locations, variables, std::vector<size_t>{1, 1});
    EXPECT_THROWS_AS(getvalues.fillGeoVaLs(state, dt_begin, dt_end, geovals),
        eckit::BadParameter);
  }

  SECTION("test fillGeoVaLs") {
    // create geovals from the locations
    std::vector<size_t> nlevs = {1};
    ufo::GeoVaLs geovals(locations, state.variables(),
        std::vector<size_t>{1, 1});

    std::stringstream os;
    os << geovals;

    getvalues.fillGeoVaLs(state, dt_begin, dt_end, geovals);

    // test fillGeoVaLs
    std::vector<double> vals(geovals.nlocs());
    geovals.get(vals, "sea_ice_area_fraction");

    double missing_value = util::missingValue(vals[0]);
    std::vector<double> testvals = {1, missing_value, 0};

    EXPECT_EQUAL(vals[0], testvals[0]);
    EXPECT_EQUAL(vals[1], testvals[1]);
    EXPECT_EQUAL(vals[2], testvals[2]);
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run( argc, argv );
}
