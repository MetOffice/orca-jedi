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
    .set("variable type", "background error variance");
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

  eckit::LocalConfiguration state_config;
  std::vector<std::string> state_variables {"sea_ice_area_fraction"};
  state_config.set("state variables", state_variables);
  state_config.set("date", "2021-06-30T00:00:00Z");
  OrcaStateParameters params;
  eckit::LocalConfiguration surf_var_conf;
  surf_var_conf.set("levels", 1);

  SECTION("test state parameters") {
    state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
    state_config.set("nemo error field file", "../Data/orca2_t_bkg_var.nc");
    params.validateAndDeserialize(state_config);
    EXPECT(params.nemoFieldFile.value() ==
        state_config.getString("nemo field file"));
    EXPECT(params.errorFieldFile.value() ==
        state_config.getString("nemo error field file"));
    EXPECT(params.analyticInit.value().value_or(true));
    auto datetime = static_cast<util::DateTime>(state_config.getString("date"));
    EXPECT(params.date.value() == datetime);
    EXPECT(params.stateVariables.value()[0].name() == state_variables[0]);
  }

  SECTION("test constructor") {
    oops::Variables oops_vars;
    for (auto && name : state_variables) {
      oops_vars.push_back(oops::Variable(name, surf_var_conf));
    }
    util::DateTime datetime("2021-06-30T00:00:00Z");
    State state(geometry, oops_vars, datetime);
  }

  SECTION("test subset copy constructor") {
    oops::Variables oops_vars({oops::Variable("sea_ice_area_fraction", surf_var_conf),
        oops::Variable("sea_surface_foundation_temperature", surf_var_conf)});
    util::DateTime datetime("2021-06-30T00:00:00Z");
    State state(geometry, oops_vars, datetime);

    State state_copy(oops::Variables({oops::Variable("sea_ice_area_fraction", surf_var_conf)}),
                     state);

    EXPECT_THROWS_AS(
        State(oops::Variables({oops::Variable("NOT_IN_SRC_STATE", surf_var_conf)}), state),
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
  state_config.set("nemo error field file", "../Data/orca2_t_bkg_var.nc");
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
  SECTION("test state to fieldset") {
    atlas::FieldSet statefset = atlas::FieldSet();
    state.State::toFieldSet(statefset);
    EXPECT(statefset[0].name() == "sea_ice_area_fraction");
  }
  SECTION("test set_gmask") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 10);
    config2.set("initialise extra fields", true);
    Geometry geometry2(config2, eckit::mpi::comm());
    params.validateAndDeserialize(state_config);
    State state(geometry2, params);
    state_config.set("set gmask", true);
    // Count number of points in gmask.
    atlas::FieldSet extraFields;
    extraFields = geometry2.extraFields();
    int gmask_sum = 0;
    for (atlas::Field field : extraFields) {
      std::string fieldname = field.name();
      std::cout << "extraFields field name: " << fieldname << std::endl;
      if (fieldname == "gmask") {
        auto field_view = atlas::array::make_view<int32_t, 2>(field);
        int field_sum = 0;
        for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
          for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
            field_sum += field_view(j, k);
          }
        }
        std::cout << fieldname << " values sum " << field_sum << std::endl;
        gmask_sum = field_sum;
      }
    }
    // The sum should be equal to the number of ocean points in the orca2 grid
    // minus ghost points and some additional masked points needed for BUMP to work.
    EXPECT(gmask_sum == 264600);
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
