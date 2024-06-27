/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include<sstream>
#include<cmath>

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

const double ATOL = 1e-6;

//-----------------------------------------------------------------------------

CASE("test basic interpolator") {
  int nlevs = 3;
  eckit::LocalConfiguration config;
  config.set("grid name", "ORCA2_T");
  config.set("number levels", nlevs);

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
  interp_conf.set("type", "unstructured-bilinear-lonlat");
  interp_conf.set("non_linear", "missing-if-all-missing-real32");
  eckit::LocalConfiguration interpolator_conf;
  interpolator_conf.set("atlas-interpolator", interp_conf);

  OrcaInterpolatorParameters params;
  params.validateAndDeserialize(interpolator_conf);

  std::vector<double> lons({0, 120, 270});
  std::vector<double> lats({88, 0, 30});
  const int nlocs = 3;

  // create a state from the test data
  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {
      "sea_ice_area_fraction",
      "sea_surface_foundation_temperature",
      "sea_water_potential_temperature"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2021-06-30T00:00:00Z");
  state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
  state_config.set("variance field file", "../Data/orca2_t_bkg_var.nc");
  OrcaStateParameters stateParams;
  stateParams.validateAndDeserialize(state_config);
  State state(geometry, stateParams);

  SECTION("test interpolator succeeds even with no locations") {
    Interpolator interpolator(interpolator_conf, geometry, {}, {});
  }

  Interpolator interpolator(interpolator_conf, geometry, lats, lons);

  SECTION("test interpolator.apply fails missing variable") {
    oops::Variables variables{{oops::Variable{"NOTAVARIABLE"}}};
    std::vector<double> vals(3);
    std::vector<bool> mask(3, true);
    EXPECT_THROWS_AS(interpolator.apply(variables, state, mask, vals),
        eckit::BadParameter);
  }

  SECTION("test interpolator.apply") {
    // two variables at n locations
    std::vector<double> vals(2*nlocs);
    std::vector<bool> mask(nlocs, true);
    interpolator.apply(oops::Variables{{oops::Variable{"sea_ice_area_fraction"},
        oops::Variable{"sea_surface_foundation_temperature"}}}, state, mask, vals);

    double missing_value = util::missingValue<double>();
    std::vector<double> testvals = {1, missing_value, 0, 18.4888916016,
                                    missing_value, 18.1592999503};

    for (size_t i=0; i < testvals.size(); ++i) {
      std::cout << "vals[" << i << "] " << std::setprecision(12) << vals[i]
                << " testvals[" << i << "] " << testvals[i] << std::endl;
      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
    }
  }
  SECTION("test interpolator.apply multiple levels") {
    std::vector<double> vals(nlevs*nlocs);
    std::vector<bool> mask(nlocs, true);
    interpolator.apply(oops::Variables{{oops::Variable{"sea_water_potential_temperature"}}},
                                       state, mask, vals);

    double missing_value = util::missingValue<double>();
    std::vector<double> testvals = {
      18.4888916016, missing_value, 18.1592999503,
      17.9419364929, missing_value, 17.75000288,
      missing_value, missing_value, missing_value};

    for (size_t i=0; i < testvals.size(); ++i) {
      std::cout << "vals[" << i << "] " << std::setprecision(12) << vals[i]
                << " testvals[" << i << "] " << testvals[i] << std::endl;
      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
    }
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
