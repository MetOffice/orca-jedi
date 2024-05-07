/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "atlas/library/Library.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "oops/util/parameters/Parameter.h"
#include "orca-jedi/state/StateParameters.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic state") {
  eckit::LocalConfiguration config;
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 2);

  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
    .set("field precision", "double")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
    .set("field precision", "double")
    .set("nemo field name", "sic_tot_var")
    .set("model space", "surface")
    .set("variable type", "background variance");
  nemo_var_mappings[2].set("name", "sea_surface_foundation_temperature")
    .set("field precision", "double")
    .set("nemo field name", "votemper")
    .set("model space", "surface");
  nemo_var_mappings[3].set("name", "sea_water_potential_temperature")
    .set("field precision", "double")
    .set("nemo field name", "votemper")
    .set("model space", "volume");
  config.set("nemo variables", nemo_var_mappings);
  Geometry geometry(config, eckit::mpi::comm());
  const std::vector<int> channels{};

  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {"sea_ice_area_fraction"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2021-06-30T00:00:00Z");
  OrcaStateParameters params;

  SECTION("test state parameters") {
    state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
    state_config.set("variance field file", "../Data/orca2_t_bkg_var.nc");
    params.validateAndDeserialize(state_config);
    EXPECT(params.nemoFieldFile.value() ==
        state_config.getString("nemo field file"));
    EXPECT(params.varianceFieldFile.value() ==
        state_config.getString("variance field file"));
    EXPECT(params.analyticInit.value().value_or(true));
    auto datetime = static_cast<util::DateTime>(state_config.getString("date"));
    EXPECT(params.date.value() == datetime);
    EXPECT(params.stateVariables.value()[0] == state_variables[0]);
  }

  SECTION("test constructor") {
    oops::Variables oops_vars(state_variables, channels);
    util::DateTime datetime("2021-06-30T00:00:00Z");
    State state(geometry, oops_vars, datetime);
  }

  SECTION("test subset copy constructor") {
    oops::Variables oops_vars({"sea_ice_area_fraction", "sea_surface_foundation_temperature"},
        channels);
    util::DateTime datetime("2021-06-30T00:00:00Z");
    State state(geometry, oops_vars, datetime);

    State state_copy(oops::Variables({"sea_ice_area_fraction"}, channels), state);

    EXPECT_THROWS_AS(State(oops::Variables({"NOT_IN_SRC_STATE"}, channels), state),
        eckit::BadValue);
  }

  SECTION("test constructor from config analytic initialisation") {
    state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
    state_config.set("analytic initialisation", true);
    params.validateAndDeserialize(state_config);
    State state(geometry, params);
    EXPECT_EQUAL(state.norm<double>("sea_ice_area_fraction"), 0);
  }

  state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
  state_config.set("variance field file", "../Data/orca2_t_bkg_var.nc");
  state_config.set("output nemo field file", "../testoutput/orca2_t_output.nc");
  params.validateAndDeserialize(state_config);
  State state(geometry, params);
  double iceNorm = 0.0032018269;
  SECTION("test constructor from state") {
    bool has_missing = state.stateFields()["sea_ice_area_fraction"].metadata()
      .has("missing_value");
    EXPECT_EQUAL(true, has_missing);
    std::cout << std::setprecision(8) << state.norm<double>("sea_ice_area_fraction")
              << std::setprecision(8) << iceNorm << std::endl;
    EXPECT(std::abs(state.norm<double>("sea_ice_area_fraction") - iceNorm) < 1e-6);
  }
  SECTION("test state read") {
    state.read(params);
    EXPECT(std::abs(state.norm<double>("sea_ice_area_fraction") - iceNorm) < 1e-6);
  }
  SECTION("test stateCopy") {
    State stateCopy(state);
    EXPECT(std::abs(stateCopy.norm<double>("sea_ice_area_fraction") - iceNorm) < 1e-6);
  }
  SECTION("test state write") {
    state.write(params);
  }

  SECTION("test state getField") {
    atlas::Field field = state.getField(0);
    std::cout << field.name() << std::endl;
    EXPECT(field.name() == "sea_ice_area_fraction");
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
