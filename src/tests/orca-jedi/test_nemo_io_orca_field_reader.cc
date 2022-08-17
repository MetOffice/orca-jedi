/*
 * (C) British Crown Copyright 2020-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/field/Field.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "eckit/exception/Exceptions.h"

#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test read_surf_var reads field array view") {
  eckit::PathName test_data_path("../Data/orca2_t_nemo.nc");

  atlas::OrcaGrid grid("ORCA2_T");
  auto meshgen_config = grid.meshgenerator();

  atlas::MeshGenerator meshgen(meshgen_config);
  auto partitioner_config = grid.partitioner();
  partitioner_config.set("type", "serial");
  auto partitioner = atlas::grid::Partitioner(partitioner_config);
  auto mesh = meshgen.generate(grid, partitioner);
  auto funcSpace = atlas::functionspace::NodeColumns(mesh);
  atlas::Field field(funcSpace.createField<double>(
                        atlas::option::name("iiceconc") |
                        atlas::option::levels(1)));

  auto field_view = atlas::array::make_view<double, 2>(field);

  NemoFieldReader field_reader(test_data_path);
  field_reader.read_surf_var("iiceconc", 0, field_view);

  atlas::Field field_2(funcSpace.createField<double>(
                        atlas::option::name("iiceconc_2") |
                        atlas::option::levels(1)));
  auto field_view_2 = atlas::array::make_view<double, 2>(field_2);

  field_reader.read_surf_var("iiceconc", mesh, 0, field_view_2);
  auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

  auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
  for (int i = 0; i<field_view.size(); ++i) {
    if (ghost(i)) continue;
    if (field_view(i, 0) != field_view_2(i, 0)) {
        std::cout << "mismatch: "
                  << " ij(" << i << ", 0) " << ij(i, 0)
                  << " ij(" << i << ", 1) " << ij(i, 1)
                  << " 1 proc " << field_view(i, 0) << " "
                  << " with mesh " << field_view_2(i, 0) << std::endl;
    }
    EXPECT_EQUAL(field_view(i, 0), field_view_2(i, 0));
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run( argc, argv );
}
