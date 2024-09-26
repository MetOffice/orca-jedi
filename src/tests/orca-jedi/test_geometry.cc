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
#include "atlas/mesh.h"

#include "orca-jedi/geometry/Geometry.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic geometry") {
  EXPECT(eckit::system::Library::exists("atlas-orca"));

  eckit::LocalConfiguration config;
  std::vector<eckit::LocalConfiguration> nemo_var_mappings(5);
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
  nemo_var_mappings[4].set("name", "depth")
    .set("nemo field name", "nav_lev")
    .set("model space", "vertical");
  config.set("nemo variables", nemo_var_mappings);
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 10);
  Geometry geometry(config, eckit::mpi::comm());

  SECTION("test geometry variable names") {
    EXPECT_THROWS_AS(geometry.nemo_var_name("NOTAVARIABLE"), eckit::BadValue);
    EXPECT(geometry.variable_in_variable_type("sea_ice_area_fraction",
                                              "background"));
    EXPECT(geometry.variable_in_variable_type("sea_ice_area_fraction_error",
                                              "background error variance"));
    EXPECT(!geometry.variable_in_variable_type("sea_ice_area_fraction_error",
                                               "background"));
    EXPECT(!geometry.variable_in_variable_type("sea_ice_area_fraction",
                                               "background error variance"));
  }

  SECTION("test geometry variable sizes") {
    oops::Variables oops_vars{{oops::Variable{"sea_ice_area_fraction"},
                               oops::Variable{"sea_water_potential_temperature"}}};
    auto varsizes = geometry.variableSizes(oops_vars);
    EXPECT_EQUAL(varsizes.size(), 2);
    EXPECT_EQUAL(varsizes[0], 1);
    EXPECT_EQUAL(varsizes[1], 10);
    oops::Variables not_vars{{oops::Variable{"NOTAVARIBLE"}}};
    EXPECT_THROWS_AS(geometry.variableSizes(not_vars), eckit::BadValue);
  }

  SECTION("test geometry variable NEMO model spaces") {
    oops::Variables oops_vars{{oops::Variable{"sea_ice_area_fraction"},
      oops::Variable{"sea_water_potential_temperature"}, oops::Variable{"depth"}}};
    auto varsizes = geometry.variableNemoSpaces(oops_vars);
    EXPECT_EQUAL(varsizes.size(), 3);
    EXPECT_EQUAL(varsizes[0], "surface");
    EXPECT_EQUAL(varsizes[1], "volume");
    EXPECT_EQUAL(varsizes[2], "vertical");

    oops::Variables not_vars{{oops::Variable{"NOTAVARIBLE"}}};
    EXPECT_THROWS_AS(geometry.variableNemoSpaces(not_vars), eckit::BadValue);

    eckit::LocalConfiguration bad_config;
    std::vector<eckit::LocalConfiguration> bad_mappings(1);
    bad_mappings[0].set("name", "sea_ice_area_fraction")
      .set("nemo field name", "iiceconc")
      .set("model space", "NONSENSICAL");
    bad_config.set("nemo variables", bad_mappings);
    bad_config.set("grid name", "ORCA2_T");
    bad_config.set("number levels", 10);
    Geometry bad_geometry(bad_config, eckit::mpi::comm());
    EXPECT_THROWS_AS(bad_geometry.variableNemoSpaces(oops_vars),
      eckit::BadValue);
  }

  SECTION("test geometry lats and lons contain no ghost points") {
    std::vector<double> lats;
    std::vector<double> lons;
    geometry.latlon(lats, lons, false);
    const auto lonlat = atlas::array::make_view<double, 2>(
        geometry.functionSpace().lonlat());
    const auto ghosts = atlas::array::make_view<int32_t, 1>(
        geometry.mesh().nodes().ghost());
    // ghost points on orca grids appear more than once in the mesh due to the
    // "seam" at the periodic boundaries. If these points are present in the
    // geometry they ought to appear more than once
    bool ghostsExist = false;
    std::vector<size_t> appearances(ghosts.size(), 0);
    for (size_t iElem = 0; iElem < ghosts.size(); ++iElem) {
      for (size_t iLoc = 0; iLoc < lats.size(); ++iLoc) {
        // skip 260, 70 as this is special case in the ORCA2_T grid where there
        // are many points overlapping
        if (lons[iLoc] - 260 < 1e-6 && lats[iLoc] - 70 < 1e-6)
          continue;
        if (lons[iLoc] == lonlat(iElem, 0) && lats[iLoc] == lonlat(iElem, 1)) {
          appearances[iElem]++;
          if (appearances[iElem] > 1) {
            std::cout << "lons[" << iLoc <<"] " << std::setprecision(12)
                      << lons[iLoc] << " lats[" << iLoc << "] "
                      << std::setprecision(12) << lats[iLoc]
                      << " ghosts(" << iElem << ") " << ghosts(iElem)
                      << " appearances[" << iElem << "] "
                      << appearances[iElem] << std::endl;
            ghostsExist = true;
          }
        }
      }
    }
    EXPECT(!ghostsExist);
  }

  SECTION("test geometry extra methods") {
    EXPECT(geometry.levelsAreTopDown());
    EXPECT(geometry.distributionType() == "serial");
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
