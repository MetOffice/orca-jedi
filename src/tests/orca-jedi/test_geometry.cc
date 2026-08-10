/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/system/LibraryManager.h"

#include "oops/base/Variables.h"

#include "atlas/array.h"
#include "atlas/library/Library.h"
#include "atlas/mesh.h"

#include "orca-jedi/geometry/Geometry.h"
#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test basic geometry") {
  EXPECT(eckit::system::LibraryManager::exists("atlas-orca"));

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

  SECTION("test geometry extrafields") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 10);
    config2.set("initialise extra fields", true);
    Geometry geometry2(config2, eckit::mpi::comm());
    atlas::FieldSet extraFields;
    extraFields = geometry2.extraFields();
    std::vector<std::string> extraFieldNames {
        "vunit",
        "owned",
        "gmask",
        "area"};
    size_t num_matches = 0;
    std::cout << "extraField list: ";
    for (atlas::Field field : extraFields) {
      std::string fieldname = field.name();
      std::cout << fieldname << " ";
      for (std::string fieldnametest : extraFieldNames) {
        if ( fieldnametest == fieldname ) {
          num_matches++;
        }
      }
    }
    std::cout << std::endl;
    std::cout << "Number of extraFields " << extraFields.size()
              << std::endl;
    std::cout << "Number matches to expected names " << num_matches
              << std::endl;
    EXPECT(static_cast<size_t>(extraFields.size()) == num_matches);
    EXPECT(extraFieldNames.size() == num_matches);
  }

  SECTION("test volume_mask not created without land sea mask config") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 10);
    config2.set("initialise extra fields", true);
    Geometry geometry2(config2, eckit::mpi::comm());
    const atlas::FieldSet& ef = geometry2.extraFields();
    EXPECT(!ef.has("volume_mask"));
  }

  SECTION("test volume_mask is independent of BUMP extra fields") {
    // No "initialise extra fields": BUMP fields (gmask etc.) are not built,
    // but set_volume_mask must still be able to create volume_mask.
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    Geometry geometry2(config2, eckit::mpi::comm());

    atlas::Field field =
        geometry2.functionSpace().createField<double>(
            atlas::option::name("test_field")
            | atlas::option::levels(3));
    const double fill = -999.0;
    field.metadata().set("missing_value", fill);
    field.metadata().set("missing_value_type", "equals");
    auto fview = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < 3; ++k) { fview(j, k) = 1.0; }
    }

    geometry2.set_volume_mask(field);

    const atlas::FieldSet& ef = geometry2.extraFields();
    EXPECT(ef.has("volume_mask"));
    EXPECT(!ef.has("gmask"));
  }

  SECTION("test set_volume_mask creates and populates volume_mask") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    config2.set("initialise extra fields", true);
    Geometry geometry2(config2, eckit::mpi::comm());

    // volume_mask should not exist yet
    EXPECT(!geometry2.extraFields().has("volume_mask"));

    // Create a field with missing values at certain (node,level) entries
    atlas::Field field =
        geometry2.functionSpace().createField<double>(
            atlas::option::name("test_field")
            | atlas::option::levels(3));
    const double fill = -999.0;
    field.metadata().set("missing_value", fill);
    field.metadata().set("missing_value_type", "equals");

    auto fview = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < 3; ++k) {
        fview(j, k) = 1.0;  // valid everywhere initially
      }
    }

    // Find a non-ghost node to mark as missing at level 1
    auto ghost = atlas::array::make_view<int32_t, 1>(
        geometry2.mesh().nodes().ghost());
    atlas::idx_t testNode = -1;
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      if (ghost(j) == 0) { testNode = j; break; }
    }
    EXPECT(testNode >= 0);
    fview(testNode, 1) = fill;  // mark level 1 as missing

    // set_volume_mask should create and populate the field
    geometry2.set_volume_mask(field);

    // volume_mask should now exist
    EXPECT(geometry2.extraFields().has("volume_mask"));
    auto vm = atlas::array::make_view<int32_t, 2>(
        geometry2.extraFields().field("volume_mask"));
    EXPECT(vm.shape(1) == 3);

    // testNode: levels 0,2 should be ocean (1), level 1 masked (0)
    EXPECT(vm(testNode, 0) == 1);
    EXPECT(vm(testNode, 1) == 0);
    EXPECT(vm(testNode, 2) == 1);
  }

  SECTION("test volume_mask is ocean on ghost nodes when unmasked") {
    // Halo-aware: with no missing values, every node - including ghost/halo
    // nodes - must be ocean (1). The old implementation force-zeroed ghost
    // and edge nodes; the hardened implementation must not.
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    Geometry geometry2(config2, eckit::mpi::comm());

    atlas::Field field =
        geometry2.functionSpace().createField<double>(
            atlas::option::name("test_field")
            | atlas::option::levels(3));
    const double fill = -999.0;
    field.metadata().set("missing_value", fill);
    field.metadata().set("missing_value_type", "equals");
    auto fview = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < 3; ++k) { fview(j, k) = 1.0; }
    }

    geometry2.set_volume_mask(field);

    auto vm = atlas::array::make_view<int32_t, 2>(
        geometry2.extraFields().field("volume_mask"));
    // Every node ocean, and at least one ghost node exists to exercise the halo.
    auto ghost = atlas::array::make_view<int32_t, 1>(
        geometry2.mesh().nodes().ghost());
    bool sawGhost = false;
    bool allOcean = true;
    for (atlas::idx_t j = 0; j < vm.shape(0); ++j) {
      if (ghost(j)) sawGhost = true;
      for (atlas::idx_t k = 0; k < vm.shape(1); ++k) {
        if (vm(j, k) != 1) allOcean = false;
      }
    }
    EXPECT(sawGhost);
    EXPECT(allOcean);
  }

  SECTION("test set_volume_mask_from_bitmask marks land by value") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    Geometry geometry2(config2, eckit::mpi::comm());

    // Explicit bitmask field: 1 = ocean, 0 = land (NEMO tmask convention).
    atlas::Field field =
        geometry2.functionSpace().createField<double>(
            atlas::option::name("tmask")
            | atlas::option::levels(3));
    auto fview = atlas::array::make_view<double, 2>(field);
    auto ghost = atlas::array::make_view<int32_t, 1>(
        geometry2.mesh().nodes().ghost());
    atlas::idx_t oceanNode = -1;
    atlas::idx_t landNode = -1;
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      const bool isGhost = ghost(j) != 0;
      // Make even owned nodes ocean, odd owned nodes land at level 0.
      for (atlas::idx_t k = 0; k < 3; ++k) { fview(j, k) = 1.0; }
      if (!isGhost && (j % 2 == 1)) {
        fview(j, 0) = 0.0;  // land at surface
        if (landNode < 0) landNode = j;
      } else if (!isGhost && oceanNode < 0) {
        oceanNode = j;
      }
    }
    EXPECT(oceanNode >= 0);
    EXPECT(landNode >= 0);

    // land value defaults to 0.
    geometry2.set_volume_mask_from_bitmask(field);

    EXPECT(geometry2.extraFields().has("volume_mask"));
    auto vm = atlas::array::make_view<int32_t, 2>(
        geometry2.extraFields().field("volume_mask"));
    EXPECT(vm(oceanNode, 0) == 1);
    EXPECT(vm(landNode, 0) == 0);
    EXPECT(vm(landNode, 1) == 1);  // only surface flagged land
  }

  SECTION("test set_volume_mask_from_bitmask honours land value") {
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    Geometry geometry2(config2, eckit::mpi::comm());

    // Inverted convention: 1 = land, 2 = ocean. land value = 1.
    atlas::Field field =
        geometry2.functionSpace().createField<float>(
            atlas::option::name("mask")
            | atlas::option::levels(3));
    auto fview = atlas::array::make_view<float, 2>(field);
    auto ghost = atlas::array::make_view<int32_t, 1>(
        geometry2.mesh().nodes().ghost());
    atlas::idx_t landNode = -1;
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < 3; ++k) { fview(j, k) = 2.0f; }  // ocean
      if (ghost(j) == 0 && landNode < 0) {
        landNode = j;
        fview(j, 2) = 1.0f;  // land at deepest level
      }
    }
    EXPECT(landNode >= 0);

    geometry2.set_volume_mask_from_bitmask(field, 1.0);

    auto vm = atlas::array::make_view<int32_t, 2>(
        geometry2.extraFields().field("volume_mask"));
    EXPECT(vm(landNode, 0) == 1);
    EXPECT(vm(landNode, 1) == 1);
    EXPECT(vm(landNode, 2) == 0);  // masked where value == land value (1)
  }

  SECTION("test volume_mask halo exchange propagates to ghost nodes") {
    // Mask every owned node at level 0. After the internal haloExchange, every
    // node (including ghost/halo nodes, which mirror owned nodes) must be
    // masked at level 0.
    eckit::LocalConfiguration config2;
    config2.set("nemo variables", nemo_var_mappings);
    config2.set("grid name", "ORCA2_T");
    config2.set("number levels", 3);
    Geometry geometry2(config2, eckit::mpi::comm());

    atlas::Field field =
        geometry2.functionSpace().createField<double>(
            atlas::option::name("test_field")
            | atlas::option::levels(3));
    const double fill = -999.0;
    field.metadata().set("missing_value", fill);
    field.metadata().set("missing_value_type", "equals");
    auto fview = atlas::array::make_view<double, 2>(field);
    auto ghost = atlas::array::make_view<int32_t, 1>(
        geometry2.mesh().nodes().ghost());
    for (atlas::idx_t j = 0; j < fview.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < 3; ++k) {
        // Missing at level 0 on owned nodes; ocean elsewhere. Ghost node
        // field values are intentionally left as-is to prove the mask on
        // ghosts comes from haloExchange, not from the ghost field values.
        fview(j, k) = (k == 0 && ghost(j) == 0) ? fill : 1.0;
      }
    }

    geometry2.set_volume_mask(field);

    auto vm = atlas::array::make_view<int32_t, 2>(
        geometry2.extraFields().field("volume_mask"));
    for (atlas::idx_t j = 0; j < vm.shape(0); ++j) {
      EXPECT(vm(j, 0) == 0);  // masked on all nodes, incl. ghosts, via halo
      EXPECT(vm(j, 1) == 1);
      EXPECT(vm(j, 2) == 1);
    }
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
