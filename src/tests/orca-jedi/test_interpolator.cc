/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include<sstream>
#include<cmath>
#include <map>

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

struct InterpTestSettingsFixture {
 public:
  eckit::LocalConfiguration geometry_config;
  eckit::LocalConfiguration state_config;
  eckit::LocalConfiguration interpolator_config;
  size_t nlocs, nlevs;
  std::vector<double> lons;
  std::vector<double> lats;
  oops::Variables surf_vars, vol_vars;
  std::vector<double> surf_values;
  std::vector<double> vol_values;
};

const double missing_value = util::missingValue<double>();

//-----------------------------------------------------------------------------

CASE("test  interpolator") {
  std::map<std::string, InterpTestSettingsFixture> settings_map{
      {"ORCA2_T", InterpTestSettingsFixture()}, {"AMM1", InterpTestSettingsFixture()}};

  // ORCA2_T settings
  {
    settings_map["ORCA2_T"].nlocs = 3;
    settings_map["ORCA2_T"].nlevs = 3;
    settings_map["ORCA2_T"].geometry_config.set("grid name", "ORCA2_T");
    settings_map["ORCA2_T"].geometry_config.set("number levels", settings_map["ORCA2_T"].nlevs);

    std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
    nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
      .set("nemo field name", "iiceconc")
      .set("model space", "surface");
    nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
      .set("nemo field name", "sic_tot_var")
      .set("model space", "surface")
      .set("variable type", "background error variance");
    nemo_var_mappings[2].set("name", "sea_surface_foundation_temperature")
      .set("nemo field name", "votemper")
      .set("model space", "surface");
    nemo_var_mappings[3].set("name", "sea_water_potential_temperature")
      .set("nemo field name", "votemper")
      .set("model space", "volume");
    settings_map["ORCA2_T"].geometry_config.set("nemo variables", nemo_var_mappings);
    eckit::LocalConfiguration interp_conf;
    interp_conf.set("type", "unstructured-bilinear-lonlat");
    interp_conf.set("non_linear", "missing-if-all-missing-real32");
    settings_map["ORCA2_T"].interpolator_config.set("atlas-interpolator", interp_conf);

    std::vector<std::string> state_variables {
        "sea_ice_area_fraction",
        "sea_surface_foundation_temperature",
        "sea_water_potential_temperature"};
    settings_map["ORCA2_T"].state_config.set("state variables", state_variables);
    settings_map["ORCA2_T"].state_config.set("date", "2021-06-30T00:00:00Z");
    settings_map["ORCA2_T"].state_config.set("nemo field file", "../Data/orca2_t_nemo.nc");
    settings_map["ORCA2_T"].state_config.set("nemo error field file", "../Data/orca2_t_bkg_var.nc");

    settings_map["ORCA2_T"].lons = std::vector<double>{0, 120, 270};
    settings_map["ORCA2_T"].lats = std::vector<double>{88, 0, 30};

    settings_map["ORCA2_T"].surf_vars = oops::Variables{{
        oops::Variable{"sea_ice_area_fraction"},
        oops::Variable{"sea_surface_foundation_temperature"}}};
    settings_map["ORCA2_T"].surf_values = std::vector<double>{
        1,             missing_value, 0,
        18.4888877869, missing_value, 18.1592998505};

    settings_map["ORCA2_T"].vol_vars = oops::Variables{
        {oops::Variable{"sea_water_potential_temperature"}}};
    settings_map["ORCA2_T"].vol_values = std::vector<double>{
        18.4888877869, missing_value, 18.1592998505,
        18           , missing_value, 17.7500019073,
        missing_value, missing_value, missing_value};
  }

  // AMM1 settings
  {
    settings_map["AMM1"].nlocs = 3;
    settings_map["AMM1"].nlevs = 3;
    settings_map["AMM1"].geometry_config.set("grid name", "../Data/amm1_atlas_grid_spec.yaml");
    settings_map["AMM1"].geometry_config.set("number levels", settings_map["AMM1"].nlevs);

    std::vector<eckit::LocalConfiguration> nemo_var_mappings(3);
    nemo_var_mappings[0].set("name", "sea_surface_height_anomaly")
      .set("nemo field name", "sossheig")
      .set("model space", "surface");
    nemo_var_mappings[1].set("name", "mass_concentration_of_chlorophyll_in_sea_water")
      .set("nemo field name", "CHL")
      .set("model space", "surface");
    nemo_var_mappings[2].set("name", "sea_water_potential_temperature")
      .set("nemo field name", "votemper")
      .set("model space", "volume");
    settings_map["AMM1"].geometry_config.set("nemo variables", nemo_var_mappings);
    eckit::LocalConfiguration interp_conf;
    interp_conf.set("type", "unstructured-bilinear-lonlat");
    interp_conf.set("non_linear", "missing-if-all-missing-real32");
    settings_map["AMM1"].interpolator_config.set("atlas-interpolator", interp_conf);

    std::vector<std::string> state_variables {
        "sea_surface_height_anomaly",
        "mass_concentration_of_chlorophyll_in_sea_water",
        "sea_water_potential_temperature"};
    settings_map["AMM1"].state_config.set("state variables", state_variables);
    settings_map["AMM1"].state_config.set("date", "2021-06-30T00:00:00Z");
    settings_map["AMM1"].state_config.set("nemo field file", "../Data/amm1_nemo.nc");
    settings_map["AMM1"].state_config.set("nemo error field file", "../Data/amm1_nemo.nc");

    settings_map["AMM1"].lons = std::vector<double>{-17.5, -6.78, -16.1};
    settings_map["AMM1"].lats = std::vector<double>{58.16, 58.91, 63.55};

    settings_map["AMM1"].surf_vars = oops::Variables{{
        oops::Variable{"sea_surface_height_anomaly"},
        oops::Variable{"mass_concentration_of_chlorophyll_in_sea_water"}}};
    settings_map["AMM1"].surf_values = std::vector<double>{
        -0.470139563084, -0.286416769028, -0.433749824762,
        0.345775008202,  0.83959633112,  2.09327220917};

    settings_map["AMM1"].vol_vars = oops::Variables{
        {oops::Variable{"sea_water_potential_temperature"}}};
    settings_map["AMM1"].vol_values = std::vector<double>{
        11.9501609802, 13.9884538651, 10.1916904449,
        11.9500904083, 13.8567619324, 10.1912517548,
        11.9499549866, 13.7289009094, 10.1908035278};
  }

  for (const auto& [key, settings] : settings_map) {
    Geometry geometry(settings.geometry_config, eckit::mpi::comm());

    OrcaInterpolatorParameters params;
    params.validateAndDeserialize(settings.interpolator_config);

    // create a state from the test data
    OrcaStateParameters stateParams;
    stateParams.validateAndDeserialize(settings.state_config);
    State state(geometry, stateParams);

    SECTION("test " + key + " interpolator succeeds even with no locations") {
      Interpolator interpolator(settings.interpolator_config, geometry, {}, {});
    }

    Interpolator interpolator(settings.interpolator_config, geometry, settings.lats, settings.lons);

    SECTION("test " + key + " interpolator.apply fails missing variable") {
      oops::Variables variables{{oops::Variable{"NOTAVARIABLE"}}};
      std::vector<double> vals(3);
      std::vector<bool> mask(3, true);
      EXPECT_THROWS_AS(interpolator.apply(variables, state, mask, vals),
          eckit::BadParameter);
    }

    SECTION("test " + key + " interpolator.apply") {
      // two variables at n locations
      std::vector<double> vals(2*settings.nlocs);
      std::vector<bool> mask(settings.nlocs, true);

      interpolator.apply(settings.surf_vars, state, mask, vals);

      for (size_t i=0; i < settings.surf_values.size(); ++i) {
        std::cout << "vals[" << i << "] " << std::setprecision(12) << vals[i]
                  << " kgo_values[" << i << "] " << settings.surf_values[i] << std::endl;
      }
      for (size_t i=0; i < settings.surf_values.size(); ++i) {
        EXPECT(std::abs(vals[i] - settings.surf_values[i]) < ATOL);
      }
    }
    SECTION("test " + key + " interpolator.apply multiple levels") {
      std::vector<double> vals(settings.nlevs*settings.nlocs);
      std::vector<bool> mask(settings.nlocs, true);

      interpolator.apply(settings.vol_vars, state, mask, vals);

      for (size_t i=0; i < settings.vol_values.size(); ++i) {
        std::cout << "vals[" << i << "] " << std::setprecision(12) << vals[i]
                  << " kgo_values[" << i << "] " << settings.vol_values[i] << std::endl;
      }
      for (size_t i=0; i < settings.vol_values.size(); ++i) {
        EXPECT(std::abs(vals[i] - settings.vol_values[i]) < ATOL);
      }
    }
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
