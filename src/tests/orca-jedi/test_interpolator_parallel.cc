/*
 * (C) British Crown Copyright 2023 Met Office
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
#include "atlas/util/function/VortexRollup.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateParameters.h"
#include "orca-jedi/interpolator/Interpolator.h"
#include "orca-jedi/interpolator/InterpolatorParameters.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

const double ATOL = 1e-6;
const double CHECKERBOARD_MED_RANK = 1;

//-----------------------------------------------------------------------------

CASE("test serial interpolator") {
  int nlevs = 3;
  eckit::LocalConfiguration config;
  config.set("grid name", "eORCA12_T");
  config.set("number levels", nlevs);
  config.set("partitioner", "serial");

  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_surface_height_anomaly")
    .set("nemo field name", "sossheig")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_surface_height_anomaly_error")
    .set("nemo field name", "ssh_tot_std")
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
  interp_conf.set("non_linear", "missing-if-all-missing");
  eckit::LocalConfiguration interpolator_conf;
  interpolator_conf.set("atlas-interpolator", interp_conf);

  OrcaInterpolatorParameters params;
  params.validateAndDeserialize(interpolator_conf);

  int nlocs = 0;
  std::vector<double> lons, lats;
  // just off baja, stick a bunch of points as these look bad in orca12 decomp
  int nlats = 10;
  int nlons = 10;
  // -114 to -94
  double lonStart = -114;
  double lonEnd = -94;
  // 10 to 15
  double latStart = 10;
  double latEnd = 15;
  nlocs = nlats*nlons;
  lons = std::vector<double>();
  lats = std::vector<double>();
  for (int iLon = 0; iLon < nlons; ++iLon) {
    double lon = lonStart + std::floor<int>((lonEnd - lonStart)*iLon/nlons);
    for (int iLat = 0; iLat < nlats; ++iLat) {
      double lat = latStart + std::floor<int>((latEnd - latStart)*iLat/nlats);
      lons.emplace_back(lon);
      lats.emplace_back(lat);
    }
  }

  auto rollup_plus = [](const double lon, const double lat) {
    return 1 + atlas::util::function::vortex_rollup(lon, lat, 0.0);
  };

  // create a state from the test data
  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {
      "sea_surface_height_anomaly",
      "sea_surface_foundation_temperature",
      "sea_water_potential_temperature"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2021-06-30T00:00:00Z");
  state_config.set("analytic initialisation", true);
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
    oops::Variables variables({"NOTAVARIABLE"});
    std::vector<double> vals(nlocs);
    std::vector<bool> mask(nlocs);
    EXPECT_THROWS_AS(interpolator.apply(variables, state, mask, vals),
        eckit::BadParameter);
  }

  SECTION("test interpolator.apply") {
    // two variables at n locations
    std::vector<double> vals(2*nlocs);
    std::vector<bool> mask(2*nlocs);
    interpolator.apply(oops::Variables({"sea_surface_height_anomaly",
        "sea_surface_foundation_temperature"}), state, mask, vals);

    double missing_value = util::missingValue<double>();
    std::vector<double> testvals(2*nlocs, 0);
    for (int i = 0; i < nlocs; ++i) {
      testvals[i] = rollup_plus(lons[i], lats[i]);
      testvals[i+nlocs] = rollup_plus(lons[i], lats[i]);
    }

    int errCount = 0;
    for (size_t i=0; i < testvals.size(); ++i) {
      if (std::abs(vals[i] - testvals[i]) > ATOL) {
      std::cout << "[" << eckit::mpi::comm().rank() << "] lons [" << i << "] " << lons[i]
                << " lats[" << i << "] " << lats[i]
                << " vals[" << i << "] " << std::setprecision(12) << vals[i]
                << " testvals[" << i << "] " << testvals[i] << std::endl;
      ++errCount;
      }
      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
    }
    std::cout << "[" << eckit::mpi::comm().rank() << "] errCount: "
              << errCount << "/" << 2*nlocs << std::endl;
  }
  //  SECTION("test interpolator.apply multiple levels") {
  //    std::vector<double> vals(nlevs*nlocs);
  //    std::vector<bool> mask(nlevs*nlocs);
  //    interpolator.apply(oops::Variables({"sea_water_potential_temperature"}),
  //                                       state, mask, vals);

  //    double missing_value = util::missingValue<double>();
  //    std::vector<double> testvals(3*nlocs,0);
  //    if (eckit::mpi::comm().rank() == CHECKERBOARD_MED_RANK) {
  //      //testvals = std::vector<double>({0.69885330268, 0.00010463659, 1.38328833133,
  //      //                                0.69885330268, 0.00010463659, 1.38328833133,
  //      //                                0.69885330268, 0.00010463659, 1.38328833133});
  //      for (int i = 0; i < nlocs; ++i) {
  //        testvals[i] = rollup_plus(lons[i], lats[i]);
  //        testvals[i+nlocs] = rollup_plus(lons[i], lats[i]);
  //        testvals[i+2*nlocs] = rollup_plus(lons[i], lats[i]);
  //      }
  //    }

  //    for (size_t i=0; i < testvals.size(); ++i) {
  //      //std::cout << "[" << eckit::mpi::comm().rank() << "] lons [" << i << "] " << lons[i]
  //      //          << " lats[" << i << "] " << lats[i]
  //      //          << " vals[" << i << "] " << std::setprecision(12) << vals[i]
  //      //          << " testvals[" << i << "] " << testvals[i] << std::endl;
  //      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
  //    }
  //  }
}

CASE("test checkerboard interpolator") {
  int nlevs = 3;
  eckit::LocalConfiguration config;
  config.set("grid name", "eORCA12_T");
  config.set("number levels", nlevs);
  config.set("source mesh halo", 0);
  config.set("partitioner", "checkerboard");

  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_surface_height_anomaly")
    .set("nemo field name", "sossheig")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_surface_height_anomaly_error")
    .set("nemo field name", "ssh_tot_std")
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
  interp_conf.set("non_linear", "missing-if-all-missing");
  eckit::LocalConfiguration interpolator_conf;
  interpolator_conf.set("atlas-interpolator", interp_conf);

  OrcaInterpolatorParameters params;
  params.validateAndDeserialize(interpolator_conf);

  int nlocs = 0;
  std::vector<double> lons, lats;
  int nlats = 10;
  int nlons = 10;
  // -114 to -94
  double lonStart = -114;
  double lonEnd = -94;
  // 10 to 15
  double latStart = 10;
  double latEnd = 15;
  nlocs = nlats*nlons;
  lons = std::vector<double>();
  lats = std::vector<double>();
  for (int iLon = 0; iLon < nlons; ++iLon) {
    double lon = lonStart + std::floor<int>((lonEnd - lonStart)*iLon/nlons);
    for (int iLat = 0; iLat < nlats; ++iLat) {
      double lat = latStart + std::floor<int>((latEnd - latStart)*iLat/nlats);
      lons.emplace_back(lon);
      lats.emplace_back(lat);
    }
  }

  auto rollup_plus = [](const double lon, const double lat) {
    return 1 + atlas::util::function::vortex_rollup(lon, lat, 0.0);
  };

  // create a state from the test data
  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {
      "sea_surface_height_anomaly",
      "sea_surface_foundation_temperature",
      "sea_water_potential_temperature"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2021-06-30T00:00:00Z");
  state_config.set("analytic initialisation", true);
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
    oops::Variables variables({"NOTAVARIABLE"});
    std::vector<double> vals(nlocs);
    std::vector<bool> mask(nlocs);
    EXPECT_THROWS_AS(interpolator.apply(variables, state, mask, vals),
        eckit::BadParameter);
  }

  SECTION("test interpolator.apply") {
    // two variables at n locations
    std::vector<double> vals(2*nlocs);
    std::vector<bool> mask(2*nlocs);
    interpolator.apply(oops::Variables({"sea_surface_height_anomaly",
        "sea_surface_foundation_temperature"}), state, mask, vals);

    double missing_value = util::missingValue<double>();
    std::vector<double> testvals(2*nlocs, 0);
    for (int i = 0; i < nlocs; ++i) {
      testvals[i] = rollup_plus(lons[i], lats[i]);
      testvals[i+nlocs] = rollup_plus(lons[i], lats[i]);
    }

    int errCount = 0;
    for (size_t i=0; i < testvals.size(); ++i) {
      if (std::abs(vals[i] - testvals[i]) > ATOL) {
      std::cout << "[" << eckit::mpi::comm().rank() << "] lons [" << i << "] " << lons[i]
                << " lats[" << i << "] " << lats[i]
                << " vals[" << i << "] " << std::setprecision(12) << vals[i]
                << " testvals[" << i << "] " << testvals[i] << std::endl;
      ++errCount;
      }
      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
    }
    std::cout << "[" << eckit::mpi::comm().rank() << "] errCount: "
              << errCount << "/" << 2*nlocs << std::endl;
  }
  //  SECTION("test interpolator.apply multiple levels") {
  //    std::vector<double> vals(nlevs*nlocs);
  //    std::vector<bool> mask(nlevs*nlocs);
  //    interpolator.apply(oops::Variables({"sea_water_potential_temperature"}),
  //                                       state, mask, vals);

  //    double missing_value = util::missingValue<double>();
  //    std::vector<double> testvals(3*nlocs,0);
  //    if (eckit::mpi::comm().rank() == CHECKERBOARD_MED_RANK) {
  //      //testvals = std::vector<double>({0.69885330268, 0.00010463659, 1.38328833133,
  //      //                                0.69885330268, 0.00010463659, 1.38328833133,
  //      //                                0.69885330268, 0.00010463659, 1.38328833133});
  //      for (int i = 0; i < nlocs; ++i) {
  //        testvals[i] = rollup_plus(lons[i], lats[i]);
  //        testvals[i+nlocs] = rollup_plus(lons[i], lats[i]);
  //        testvals[i+2*nlocs] = rollup_plus(lons[i], lats[i]);
  //      }
  //    }

  //    for (size_t i=0; i < testvals.size(); ++i) {
  //      //std::cout << "[" << eckit::mpi::comm().rank() << "] lons [" << i << "] " << lons[i]
  //      //          << " lats[" << i << "] " << lats[i]
  //      //          << " vals[" << i << "] " << std::setprecision(12) << vals[i]
  //      //          << " testvals[" << i << "] " << testvals[i] << std::endl;
  //      EXPECT(std::abs(vals[i] - testvals[i]) < ATOL);
  //    }
  //  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
