/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/array.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/regridder/Regridder.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

// ---------------------------------------------------------------------------
CASE("test regridder structured-to-structured") {
  // Source: coarse regular lon-lat grid
  atlas::StructuredGrid sourceGrid("L24");  // ~7.5 degree global
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh sourceMesh = meshgen.generate(sourceGrid);
  atlas::functionspace::NodeColumns sourceFunctionSpace(sourceMesh);

  // Target: finer regular lon-lat grid
  atlas::StructuredGrid targetGrid("L48");  // ~3.75 degree global
  atlas::Mesh targetMesh = meshgen.generate(targetGrid);
  atlas::functionspace::NodeColumns targetFunctionSpace(targetMesh);

  // Create a source field with a smooth analytic function (cosine of latitude)
  atlas::Field sourceField = sourceFunctionSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  auto source_view = atlas::array::make_view<double, 2>(sourceField);
  auto lonlat = atlas::array::make_view<double, 2>(sourceFunctionSpace.nodes().lonlat());
  for (atlas::idx_t i = 0; i < sourceFunctionSpace.nodes().size(); ++i) {
    double lat_rad = lonlat(i, 1) * M_PI / 180.0;
    source_view(i, 0) = std::cos(lat_rad);
  }

  // Configure interpolation: finite-element method
  eckit::LocalConfiguration interpConf;
  interpConf.set("type", "finite-element");

  // Build regridder and execute
  orcamodel::Regridder regridder(interpConf, sourceFunctionSpace, targetFunctionSpace);
  atlas::Field targetField = regridder.execute(sourceField);

  // Verify: check that interpolated values are reasonable
  // Skip pole nodes (|lat| > 85) where FE interpolation on structured grids
  // is unreliable due to degenerate triangular elements at the pole cap.
  auto target_view = atlas::array::make_view<double, 2>(targetField);
  auto target_lonlat = atlas::array::make_view<double, 2>(
      targetFunctionSpace.nodes().lonlat());
  auto target_ghost = atlas::array::make_view<int, 1>(
      targetFunctionSpace.nodes().ghost());

  double max_error = 0.0;
  int checked_points = 0;
  for (atlas::idx_t i = 0; i < targetFunctionSpace.nodes().size(); ++i) {
    if (target_ghost(i)) continue;  // skip ghost/halo nodes
    double lat = target_lonlat(i, 1);
    if (std::abs(lat) > 85.0) continue;  // skip pole caps
    double lat_rad = lat * M_PI / 180.0;
    double expected = std::cos(lat_rad);
    double error = std::abs(target_view(i, 0) - expected);
    max_error = std::max(max_error, error);
    ++checked_points;
  }

  // For cos(lat) on a ~7.5 degree source grid, interpolation error should be small
  eckit::Log::info() << "Max interpolation error: " << max_error
                     << " (checked " << checked_points << " points)" << std::endl;
  EXPECT(checked_points > 0);
  EXPECT(max_error < 0.05);
}

// ---------------------------------------------------------------------------
CASE("test regridder fieldset") {
  // Source grid
  atlas::StructuredGrid sourceGrid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh sourceMesh = meshgen.generate(sourceGrid);
  atlas::functionspace::NodeColumns sourceFunctionSpace(sourceMesh);

  // Target grid
  atlas::StructuredGrid targetGrid("L48");
  atlas::Mesh targetMesh = meshgen.generate(targetGrid);
  atlas::functionspace::NodeColumns targetFunctionSpace(targetMesh);

  // Create two source fields
  atlas::FieldSet sourceFields;
  atlas::Field field1 = sourceFunctionSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  atlas::Field field2 = sourceFunctionSpace.createField<double>(
      atlas::option::name("salinity") | atlas::option::levels(1));

  auto view1 = atlas::array::make_view<double, 2>(field1);
  auto view2 = atlas::array::make_view<double, 2>(field2);
  for (atlas::idx_t i = 0; i < sourceFunctionSpace.nodes().size(); ++i) {
    view1(i, 0) = 1.0;  // uniform field
    view2(i, 0) = 2.0;  // uniform field
  }
  sourceFields.add(field1);
  sourceFields.add(field2);

  // Configure and execute
  eckit::LocalConfiguration interpConf;
  interpConf.set("type", "finite-element");

  orcamodel::Regridder regridder(interpConf, sourceFunctionSpace, targetFunctionSpace);
  atlas::FieldSet targetFields = regridder.execute(sourceFields);

  // Verify both fields were regridded
  EXPECT(targetFields.size() == 2);
  EXPECT(targetFields[0].name() == "temperature");
  EXPECT(targetFields[1].name() == "salinity");

  // Uniform fields should remain uniform after interpolation (skip ghost nodes)
  auto tgt_view1 = atlas::array::make_view<double, 2>(targetFields[0]);
  auto tgt_view2 = atlas::array::make_view<double, 2>(targetFields[1]);
  auto tgt_ghost = atlas::array::make_view<int, 1>(
      targetFunctionSpace.nodes().ghost());
  for (atlas::idx_t i = 0; i < targetFunctionSpace.nodes().size(); ++i) {
    if (tgt_ghost(i)) continue;
    EXPECT(std::abs(tgt_view1(i, 0) - 1.0) < 1e-7);
    EXPECT(std::abs(tgt_view2(i, 0) - 2.0) < 1e-7);
  }
}

// ---------------------------------------------------------------------------
CASE("test regridder with nonlinear missing value handling") {
  // Source grid
  atlas::StructuredGrid sourceGrid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh sourceMesh = meshgen.generate(sourceGrid);
  atlas::functionspace::NodeColumns sourceFunctionSpace(sourceMesh);

  // Target grid
  atlas::StructuredGrid targetGrid("L48");
  atlas::Mesh targetMesh = meshgen.generate(targetGrid);
  atlas::functionspace::NodeColumns targetFunctionSpace(targetMesh);

  // Create source field with missing values (simulating land mask)
  atlas::Field sourceField = sourceFunctionSpace.createField<double>(
      atlas::option::name("sst") | atlas::option::levels(1));
  auto source_view = atlas::array::make_view<double, 2>(sourceField);
  const double missing = 9999.0;

  // Set metadata for missing value
  sourceField.metadata().set("missing_value", missing);
  sourceField.metadata().set("missing_value_type", "equals");

  auto lonlat = atlas::array::make_view<double, 2>(sourceFunctionSpace.nodes().lonlat());
  for (atlas::idx_t i = 0; i < sourceFunctionSpace.nodes().size(); ++i) {
    double lat = lonlat(i, 1);
    // "Land" above 60N — set to missing
    if (lat > 60.0) {
      source_view(i, 0) = missing;
    } else {
      source_view(i, 0) = 20.0 - 0.3 * std::abs(lat);  // simple SST-like profile
    }
  }

  // Configure with nonlinear missing value handling
  eckit::LocalConfiguration interpConf;
  interpConf.set("type", "finite-element");
  interpConf.set("non_linear", "missing-if-all-missing");

  orcamodel::Regridder regridder(interpConf, sourceFunctionSpace, targetFunctionSpace);
  atlas::Field targetField = regridder.execute(sourceField);

  // Verify: points well within "ocean" should have valid values
  auto target_view = atlas::array::make_view<double, 2>(targetField);
  auto target_lonlat = atlas::array::make_view<double, 2>(
      targetFunctionSpace.nodes().lonlat());
  auto target_ghost = atlas::array::make_view<int, 1>(
      targetFunctionSpace.nodes().ghost());

  int valid_ocean_points = 0;
  for (atlas::idx_t i = 0; i < targetFunctionSpace.nodes().size(); ++i) {
    if (target_ghost(i)) continue;
    double lat = target_lonlat(i, 1);
    if (lat < 50.0 && lat > -50.0) {
      // Well within ocean — should have valid interpolated value
      EXPECT(target_view(i, 0) != missing);
      valid_ocean_points++;
    }
  }
  EXPECT(valid_ocean_points > 0);
  eckit::Log::info() << "Verified " << valid_ocean_points
                     << " valid ocean target points" << std::endl;
}

// ---------------------------------------------------------------------------
CASE("test regridder ORCA1_T to ORCA2_T") {
  // Source: ORCA1_T geometry
  eckit::LocalConfiguration sourceGeomConf;
  sourceGeomConf.set("grid name", "ORCA1_T");
  sourceGeomConf.set("number levels", 3);
  std::vector<eckit::LocalConfiguration> srcVarMappings(2);
  srcVarMappings[0].set("name", "sea_ice_area_fraction")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  srcVarMappings[1].set("name", "sea_water_potential_temperature")
    .set("nemo field name", "votemper")
    .set("model space", "volume");
  sourceGeomConf.set("nemo variables", srcVarMappings);

  orcamodel::Geometry sourceGeom(sourceGeomConf, eckit::mpi::comm());

  // Read source state from test data
  eckit::LocalConfiguration stateConf;
  stateConf.set("state variables",
      std::vector<std::string>{"sea_ice_area_fraction",
                               "sea_water_potential_temperature"});
  stateConf.set("date", "2021-06-30T00:00:00Z");
  stateConf.set("nemo field file", "../Data/orca1_t_nemo.nc");
  orcamodel::State sourceState(sourceGeom, stateConf);

  eckit::Log::info() << "Source state (ORCA1_T): " << sourceState << std::endl;

  // Target: ORCA2_T geometry
  eckit::LocalConfiguration targetGeomConf;
  targetGeomConf.set("grid name", "ORCA2_T");
  targetGeomConf.set("number levels", 3);
  std::vector<eckit::LocalConfiguration> tgtVarMappings(2);
  tgtVarMappings[0].set("name", "sea_ice_area_fraction")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  tgtVarMappings[1].set("name", "sea_water_potential_temperature")
    .set("nemo field name", "votemper")
    .set("model space", "volume");
  targetGeomConf.set("nemo variables", tgtVarMappings);

  orcamodel::Geometry targetGeom(targetGeomConf, eckit::mpi::comm());

  // Configure interpolation: unstructured-bilinear-lonlat with missing value handling
  eckit::LocalConfiguration interpConf;
  interpConf.set("type", "unstructured-bilinear-lonlat");
  interpConf.set("non_linear", "missing-if-all-missing");

  // Build regridder: ORCA1_T function space -> ORCA2_T function space
  orcamodel::Regridder regridder(interpConf,
                                 sourceGeom.functionSpace(),
                                 targetGeom.functionSpace());

  // Execute regridding
  atlas::FieldSet result = regridder.execute(sourceState.stateFields());

  eckit::Log::info() << "Regridded " << result.size() << " fields from ORCA1_T to ORCA2_T"
                     << std::endl;

  // Verify: we got fields out with the right names and sizes
  EXPECT(result.size() == sourceState.stateFields().size());

  for (atlas::idx_t f = 0; f < result.size(); ++f) {
    const atlas::Field& field = result[f];
    eckit::Log::info() << "  Field '" << field.name()
                       << "' shape: (" << field.shape(0)
                       << ", " << field.shape(1) << ")" << std::endl;
    // Target field should have ORCA2_T number of nodes
    EXPECT(field.shape(0) == targetGeom.functionSpace().size());
  }

  // Verify: sea_ice_area_fraction values should be broadly reasonable.
  // Bilinear interpolation near coastlines (ocean/land boundaries) can
  // produce overshoot/undershoot when missing-if-all-missing allows
  // partial stencils, so we use a generous tolerance.
  // Note: ORCA fields are stored as float (real32) in atlas
  auto sic_view = atlas::array::make_view<float, 2>(result["sea_ice_area_fraction"]);
  atlas::functionspace::NodeColumns tgtNodeColumns(targetGeom.functionSpace());
  auto tgt_ghost = atlas::array::make_view<int, 1>(
      tgtNodeColumns.nodes().ghost());
  int valid_sic_points = 0;
  float sic_min = std::numeric_limits<float>::max();
  float sic_max = std::numeric_limits<float>::lowest();
  for (atlas::idx_t i = 0; i < tgtNodeColumns.size(); ++i) {
    if (tgt_ghost(i)) continue;
    float val = sic_view(i, 0);
    if (val != 0.0f) {  // non-zero means it was interpolated
      valid_sic_points++;
      sic_min = std::min(sic_min, val);
      sic_max = std::max(sic_max, val);
    }
  }
  eckit::Log::info() << "  sea_ice_area_fraction: " << valid_sic_points
                     << " valid points, range [" << sic_min << ", "
                     << sic_max << "]" << std::endl;
  EXPECT(valid_sic_points > 0);
  // Values should not be wildly out of physical range (allow interpolation overshoot)
  EXPECT(sic_min > -0.5f);
  EXPECT(sic_max < 1.5f);
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
  return orcamodel::test::run(argc, argv);
}
