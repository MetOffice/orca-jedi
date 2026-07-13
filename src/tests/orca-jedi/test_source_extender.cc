/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <cmath>
#include <limits>
#include <vector>

#include "eckit/testing/Test.h"

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "orca-jedi/regridder/SourceExtender.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

// ---------------------------------------------------------------------------
CASE("test source extender flood fill cell-based") {
  // Create a structured mesh with NodeColumns function space
  atlas::StructuredGrid grid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  atlas::functionspace::NodeColumns funcSpace(mesh);

  // Create a source field: "ocean" where |lat| < 60, "land" elsewhere
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  const double missing = 9999.0;
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");

  auto view = atlas::array::make_view<double, 2>(field);
  auto lonlat = atlas::array::make_view<double, 2>(funcSpace.nodes().lonlat());
  auto ghost = atlas::array::make_view<int, 1>(funcSpace.nodes().ghost());

  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    double lat = lonlat(i, 1);
    if (std::abs(lat) > 60.0) {
      view(i, 0) = missing;  // land
    } else {
      view(i, 0) = 20.0 - 0.3 * std::abs(lat);  // SST-like
    }
  }

  // Count initial missing points (non-ghost)
  int initialMissing = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (view(i, 0) == missing) ++initialMissing;
  }
  eckit::Log::info() << "Initial missing (non-ghost): " << initialMissing
                     << std::endl;
  EXPECT(initialMissing > 0);

  // Apply cell-based flood fill with 2 iterations
  orcamodel::SourceExtender extender(
      mesh, funcSpace,
      /*nFloodIterations=*/2,
      /*nSmoothIterations=*/10,
      /*smoothWeightSelf=*/0.35,
      orcamodel::AdjacencyType::CellBased);
  extender.extend(field);

  // Count remaining missing points (should be fewer)
  int remainingMissing = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (view(i, 0) == missing) ++remainingMissing;
  }
  eckit::Log::info() << "Remaining missing (non-ghost): " << remainingMissing
                     << std::endl;
  EXPECT(remainingMissing < initialMissing);

  // Points that were flooded should have reasonable values (not missing, not NaN)
  int floodedValid = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    double val = view(i, 0);
    if (val != missing && std::abs(lonlat(i, 1)) > 60.0) {
      // Was originally land, now has a value
      EXPECT(!std::isnan(val));
      EXPECT(std::isfinite(val));
      ++floodedValid;
    }
  }
  eckit::Log::info() << "Flooded valid points: " << floodedValid << std::endl;
  EXPECT(floodedValid > 0);
}

// ---------------------------------------------------------------------------
CASE("test source extender flood fill edge-based") {
  atlas::StructuredGrid grid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  atlas::functionspace::NodeColumns funcSpace(mesh);

  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  const double missing = 9999.0;
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");

  auto view = atlas::array::make_view<double, 2>(field);
  auto lonlat = atlas::array::make_view<double, 2>(funcSpace.nodes().lonlat());
  auto ghost = atlas::array::make_view<int, 1>(funcSpace.nodes().ghost());

  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    double lat = lonlat(i, 1);
    if (std::abs(lat) > 60.0) {
      view(i, 0) = missing;
    } else {
      view(i, 0) = 20.0 - 0.3 * std::abs(lat);
    }
  }

  int initialMissing = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (view(i, 0) == missing) ++initialMissing;
  }

  // Apply edge-based flood fill
  orcamodel::SourceExtender extender(
      mesh, funcSpace,
      /*nFloodIterations=*/2,
      /*nSmoothIterations=*/10,
      /*smoothWeightSelf=*/0.35,
      orcamodel::AdjacencyType::EdgeBased);
  extender.extend(field);

  int remainingMissing = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (view(i, 0) == missing) ++remainingMissing;
  }
  eckit::Log::info() << "Edge-based remaining missing (non-ghost): "
                     << remainingMissing << std::endl;
  EXPECT(remainingMissing < initialMissing);
}

// ---------------------------------------------------------------------------
CASE("test source extender no-op when no missing values") {
  atlas::StructuredGrid grid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  atlas::functionspace::NodeColumns funcSpace(mesh);

  // Field with no missing values — should be unchanged
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  const double missing = 9999.0;
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");

  auto view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    view(i, 0) = 15.0;  // all valid
  }

  orcamodel::SourceExtender extender(mesh, funcSpace);
  extender.extend(field);

  // All values should remain 15.0
  auto ghost = atlas::array::make_view<int, 1>(funcSpace.nodes().ghost());
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    EXPECT(std::abs(view(i, 0) - 15.0) < 1e-10);
  }
}

// ---------------------------------------------------------------------------
CASE("test source extender skips field without missing value metadata") {
  atlas::StructuredGrid grid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  atlas::functionspace::NodeColumns funcSpace(mesh);

  // Field without missing_value metadata — should be skipped
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(1));
  auto view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    view(i, 0) = 42.0;
  }

  orcamodel::SourceExtender extender(mesh, funcSpace);
  extender.extend(field);  // should not crash

  // Value unchanged
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    EXPECT(std::abs(view(i, 0) - 42.0) < 1e-10);
  }
}

// ---------------------------------------------------------------------------
CASE("test source extender fieldset") {
  atlas::StructuredGrid grid("L24");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  atlas::functionspace::NodeColumns funcSpace(mesh);

  const double missing = 9999.0;

  // Two fields with different land patterns
  atlas::Field field1 = funcSpace.createField<double>(
      atlas::option::name("sst") | atlas::option::levels(1));
  field1.metadata().set("missing_value", missing);
  field1.metadata().set("missing_value_type", "equals");

  atlas::Field field2 = funcSpace.createField<double>(
      atlas::option::name("salinity") | atlas::option::levels(1));
  field2.metadata().set("missing_value", missing);
  field2.metadata().set("missing_value_type", "equals");

  auto view1 = atlas::array::make_view<double, 2>(field1);
  auto view2 = atlas::array::make_view<double, 2>(field2);
  auto lonlat = atlas::array::make_view<double, 2>(funcSpace.nodes().lonlat());

  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    double lat = lonlat(i, 1);
    if (std::abs(lat) > 60.0) {
      view1(i, 0) = missing;
      view2(i, 0) = missing;
    } else {
      view1(i, 0) = 20.0;
      view2(i, 0) = 35.0;
    }
  }

  atlas::FieldSet fields;
  fields.add(field1);
  fields.add(field2);

  orcamodel::SourceExtender extender(mesh, funcSpace);
  extender.extend(fields);

  // Both fields should have fewer missing values
  auto ghost = atlas::array::make_view<int, 1>(funcSpace.nodes().ghost());
  int missing1 = 0, missing2 = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (view1(i, 0) == missing) ++missing1;
    if (view2(i, 0) == missing) ++missing2;
  }
  eckit::Log::info() << "FieldSet remaining missing: sst=" << missing1
                     << " salinity=" << missing2 << std::endl;

  // Count original missing for comparison
  int origMissing = 0;
  for (atlas::idx_t i = 0; i < funcSpace.nodes().size(); ++i) {
    if (ghost(i)) continue;
    if (std::abs(lonlat(i, 1)) > 60.0) ++origMissing;
  }
  EXPECT(missing1 < origMissing);
  EXPECT(missing2 < origMissing);
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
  return orcamodel::test::run(argc, argv);
}
