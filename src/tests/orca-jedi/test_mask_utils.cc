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

#include "orca-jedi/utilities/MaskUtils.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

// Helper: create a simple structured mesh and NodeColumns function space
static atlas::functionspace::NodeColumns makeNodeColumns() {
  atlas::StructuredGrid grid("L8");
  atlas::MeshGenerator meshgen("structured");
  atlas::Mesh mesh = meshgen.generate(grid);
  return atlas::functionspace::NodeColumns(mesh);
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields surface mask broadcasts to all levels") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();
  const atlas::idx_t nLevels = 3;
  const double missing = -999.0;

  // Create a surface mask (nNodes, 1): mask out every other node
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    mask_view(j, 0) = (j % 2 == 0) ? 1 : 0;  // even=ocean, odd=land
  }

  // Create a volumetric field with 3 levels, all values = 42.0
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("temperature") | atlas::option::levels(nLevels));
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");
  auto field_view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    for (atlas::idx_t k = 0; k < nLevels; ++k) {
      field_view(j, k) = 42.0;
    }
  }

  atlas::FieldSet fields;
  fields.add(field);
  orcamodel::applyMaskToFields(mask, fields);

  // Verify: masked (odd) nodes should have missing at ALL levels
  //         valid (even) nodes should remain 42.0
  auto result = atlas::array::make_view<double, 2>(fields[0]);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    for (atlas::idx_t k = 0; k < nLevels; ++k) {
      if (j % 2 == 0) {
        EXPECT(result(j, k) == 42.0);
      } else {
        EXPECT(result(j, k) == missing);
      }
    }
  }
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields volumetric mask per-level") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();
  const atlas::idx_t nLevels = 3;
  const double missing = -999.0;

  // Create a volumetric mask (nNodes, 3):
  // level 0: all ocean (1), level 1: mask odd nodes, level 2: all masked (0)
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("volume_mask") | atlas::option::levels(nLevels));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    mask_view(j, 0) = 1;                         // level 0: all ocean
    mask_view(j, 1) = (j % 2 == 0) ? 1 : 0;     // level 1: even=ocean
    mask_view(j, 2) = 0;                         // level 2: all land
  }

  // Create field, all values = 10.0
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("salinity") | atlas::option::levels(nLevels));
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");
  auto field_view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    for (atlas::idx_t k = 0; k < nLevels; ++k) {
      field_view(j, k) = 10.0;
    }
  }

  atlas::FieldSet fields;
  fields.add(field);
  orcamodel::applyMaskToFields(mask, fields);

  auto result = atlas::array::make_view<double, 2>(fields[0]);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    // Level 0: all should remain valid
    EXPECT(result(j, 0) == 10.0);
    // Level 1: even nodes valid, odd nodes masked
    if (j % 2 == 0) {
      EXPECT(result(j, 1) == 10.0);
    } else {
      EXPECT(result(j, 1) == missing);
    }
    // Level 2: all should be masked
    EXPECT(result(j, 2) == missing);
  }
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields with float (real32) fields") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();
  const float missing = -999.0f;

  // Surface mask: mask node 0 only
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    mask_view(j, 0) = (j == 0) ? 0 : 1;
  }

  // Create a float field
  atlas::Field field = funcSpace.createField<float>(
      atlas::option::name("sst") | atlas::option::levels(1));
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");
  auto field_view = atlas::array::make_view<float, 2>(field);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    field_view(j, 0) = 25.0f;
  }

  atlas::FieldSet fields;
  fields.add(field);
  orcamodel::applyMaskToFields(mask, fields);

  auto result = atlas::array::make_view<float, 2>(fields[0]);
  EXPECT(result(0, 0) == missing);
  for (atlas::idx_t j = 1; j < nNodes; ++j) {
    EXPECT(result(j, 0) == 25.0f);
  }
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields skips fields without missing_value") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();

  // Mask that would mask everything
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    mask_view(j, 0) = 0;
  }

  // Field with no missing_value metadata
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("no_mv") | atlas::option::levels(1));
  auto field_view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    field_view(j, 0) = 7.0;
  }

  atlas::FieldSet fields;
  fields.add(field);
  orcamodel::applyMaskToFields(mask, fields);

  // Field should be unchanged (all 7.0) since it has no missing_value
  auto result = atlas::array::make_view<double, 2>(fields[0]);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    EXPECT(result(j, 0) == 7.0);
  }
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields vol mask with fewer field levels") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();
  const double missing = -999.0;

  // Volumetric mask with 5 levels
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("volume_mask") | atlas::option::levels(5));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    for (atlas::idx_t k = 0; k < 5; ++k) {
      // Mask level 2 for all nodes
      mask_view(j, k) = (k == 2) ? 0 : 1;
    }
  }

  // Field with only 2 levels (fewer than mask)
  atlas::Field field = funcSpace.createField<double>(
      atlas::option::name("shallow") | atlas::option::levels(2));
  field.metadata().set("missing_value", missing);
  field.metadata().set("missing_value_type", "equals");
  auto field_view = atlas::array::make_view<double, 2>(field);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    field_view(j, 0) = 1.0;
    field_view(j, 1) = 2.0;
  }

  atlas::FieldSet fields;
  fields.add(field);
  orcamodel::applyMaskToFields(mask, fields);

  // Field levels 0,1 are valid in the mask — should be unchanged
  auto result = atlas::array::make_view<double, 2>(fields[0]);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    EXPECT(result(j, 0) == 1.0);
    EXPECT(result(j, 1) == 2.0);
  }
}

// ---------------------------------------------------------------------------
CASE("test applyMaskToFields multiple fields mixed types") {
  auto funcSpace = makeNodeColumns();
  const atlas::idx_t nNodes = funcSpace.size();
  const double missing_d = -999.0;
  const float missing_f = -999.0f;

  // Surface mask: mask node 0
  atlas::Field mask = funcSpace.createField<int32_t>(
      atlas::option::name("mask") | atlas::option::levels(1));
  auto mask_view = atlas::array::make_view<int32_t, 2>(mask);
  for (atlas::idx_t j = 0; j < nNodes; ++j) {
    mask_view(j, 0) = (j == 0) ? 0 : 1;
  }

  // Double field
  atlas::Field field_d = funcSpace.createField<double>(
      atlas::option::name("temp") | atlas::option::levels(1));
  field_d.metadata().set("missing_value", missing_d);
  field_d.metadata().set("missing_value_type", "equals");
  auto view_d = atlas::array::make_view<double, 2>(field_d);
  for (atlas::idx_t j = 0; j < nNodes; ++j) view_d(j, 0) = 5.0;

  // Float field
  atlas::Field field_f = funcSpace.createField<float>(
      atlas::option::name("salt") | atlas::option::levels(1));
  field_f.metadata().set("missing_value", missing_f);
  field_f.metadata().set("missing_value_type", "equals");
  auto view_f = atlas::array::make_view<float, 2>(field_f);
  for (atlas::idx_t j = 0; j < nNodes; ++j) view_f(j, 0) = 35.0f;

  atlas::FieldSet fields;
  fields.add(field_d);
  fields.add(field_f);
  orcamodel::applyMaskToFields(mask, fields);

  // Node 0 should be masked in both fields
  auto res_d = atlas::array::make_view<double, 2>(fields[0]);
  auto res_f = atlas::array::make_view<float, 2>(fields[1]);
  EXPECT(res_d(0, 0) == missing_d);
  EXPECT(res_f(0, 0) == missing_f);
  // Node 1 should be valid in both
  EXPECT(res_d(1, 0) == 5.0);
  EXPECT(res_f(1, 0) == 35.0f);
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
  return orcamodel::test::run(argc, argv);
}
