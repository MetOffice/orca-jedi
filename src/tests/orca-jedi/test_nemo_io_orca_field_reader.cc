/*
 * (C) British Crown Copyright 2023 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/testing/Test.h"

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

CASE("test parallel serially distributed reads field array view") {
  eckit::PathName test_data_path("../Data/orca2_t_nemo.nc");

  atlas::OrcaGrid grid("ORCA2_T");
  auto meshgen_config = grid.meshgenerator();

  atlas::MeshGenerator meshgen(meshgen_config);
  auto partitioner_config = grid.partitioner();
  partitioner_config.set("type", "serial");
  auto partitioner = atlas::grid::Partitioner(partitioner_config);
  auto mesh = meshgen.generate(grid, partitioner);
  auto funcSpace = atlas::functionspace::NodeColumns(mesh);

  NemoFieldReader field_reader(test_data_path);
  std::vector<double> raw_data = field_reader.read_var_slice("iiceconc", 0, 0);

  atlas::Field field(funcSpace.createField<double>(
                        atlas::option::name("iiceconc") |
                        atlas::option::levels(1)));
  auto field_view = atlas::array::make_view<double, 2>(field);

  field_reader.read_surf_var("iiceconc", mesh, 0, field_view);
  auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

  auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
  for (size_t i = 0; i < field_view.size(); ++i) {
    if (ghost(i)) continue;
    if (raw_data[i] != field_view(i, 0)) {
        std::cout << "mismatch: "
                  << " ij(" << i << ", 0) " << ij(i, 0)
                  << " ij(" << i << ", 1) " << ij(i, 1)
                  << " 1 proc " << raw_data[i] << " "
                  << " with mesh " << field_view(i, 0) << std::endl;
    }
    EXPECT_EQUAL(raw_data[i], field_view(i, 0));
  }

  SECTION("volume field") {
    atlas::Field field(funcSpace.createField<double>(
                          atlas::option::name("votemper") |
                          atlas::option::levels(3)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    field_reader.read_volume_var("votemper", mesh, 0, field_view);
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    std::vector<double> raw_data;
    for (int k =0; k <3; k++) {
      raw_data = field_reader.read_var_slice("votemper", 0, k);
      for (int i = 0; i < field_view.shape(0); ++i) {
        if (ghost(i)) continue;
        if (raw_data[i] != field_view(i, k)) {
            std::cout << "mismatch: "
                      << " ij(" << i << ", 0) " << ij(i, 0)
                      << " ij(" << i << ", 1) " << ij(i, 1)
                      << " 1 proc " << raw_data[i] << " "
                      << " with mesh " << field_view(i, k) << std::endl;
        }
        EXPECT_EQUAL(raw_data[i], field_view(i, k));
      }
    }
  }

  SECTION("depth field") {
    atlas::Field field(funcSpace.createField<double>(
                          atlas::option::name("depth") |
                          atlas::option::levels(3)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    field_reader.read_vertical_var("nav_lev", mesh, field_view);
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    std::vector<double> levels{0, 10, 100};
    for (int k =0; k <3; k++) {
      for (int i = 0; i < field_view.shape(0); ++i) {
        if (ghost(i)) continue;
        if (levels[k] != field_view(i, k)) {
            std::cout << "mismatch: "
                      << " ij(" << i << ", 0) " << ij(i, 0)
                      << " ij(" << i << ", 1) " << ij(i, 1)
                      << " 1 proc " << levels[k] << " "
                      << " with mesh " << field_view(i, k) << std::endl;
        }
        EXPECT_EQUAL(levels[k], field_view(i, k));
      }
    }
  }
}

// Disable this section until jopa-bundle uses a new version of atlas
/*
CASE("test parallel domain distributed read_surf_var reads field array view") {
  // NOTE: At this time, atlas-orca is only capable of domain distribution for
  //       ORCA1 and higher resolution
  eckit::PathName test_data_path("../Data/orca1_t_nemo.nc");

  atlas::OrcaGrid grid("ORCA1_T");
  auto meshgen_config = grid.meshgenerator();

  atlas::MeshGenerator meshgen(meshgen_config);
  auto partitioner_config = grid.partitioner();
  partitioner_config.set("type", "checkerboard");
  auto partitioner = atlas::grid::Partitioner(partitioner_config);
  auto mesh = meshgen.generate(grid, partitioner);
  auto funcSpace = atlas::functionspace::NodeColumns(mesh);
  atlas::Field field(funcSpace.createField<double>(
                        atlas::option::name("iiceconc") |
                        atlas::option::levels(1)));

  auto field_view = atlas::array::make_view<double, 2>(field);

  NemoFieldReader field_reader(test_data_path);
  field_reader.read_surf_var("iiceconc", mesh, 0, field_view);

  std::vector<double> raw_data = field_reader.read_var_slice("iiceconc", 0, 0);

  int nx_halo_WE = grid.nx() + grid.haloEast() + grid.haloWest();
  int ny_halo_NS = grid.ny() + grid.haloNorth() + grid.haloSouth();

  int iy_glb_min = -grid.haloSouth();
  int ix_glb_min = -grid.haloWest();
  int glbarray_offset  = -( nx_halo_WE * iy_glb_min ) - ix_glb_min;
  int glbarray_jstride = nx_halo_WE;

  auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
  auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
  for (size_t i = 0; i<mesh.nodes().size(); ++i) {
    if (ghost(i)) continue;
    int orca_i = ij(i, 0);
    int orca_j = ij(i, 1);
    size_t raw_idx = glbarray_offset + orca_j * glbarray_jstride + orca_i;
    if (field_view(i, 0) != raw_data[raw_idx]) {
        std::cout << "mismatch: "
                  << " ij(" << i << ", 0) " << ij(i, 0)
                  << " ij(" << i << ", 1) " << ij(i, 1)
                  << " checkerboard " << field_view(i, 0) << " "
                  << " raw " << raw_data[raw_idx] << std::endl;
    }
    EXPECT_EQUAL(field_view(i, 0), raw_data[raw_idx]);
  }

  SECTION("volume field") {
    atlas::Field field(funcSpace.createField<double>(
                          atlas::option::name("votemper") |
                          atlas::option::levels(3)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    field_reader.read_volume_var("votemper", mesh, 0, field_view);
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    std::vector<double> raw_data;
    for (int k =0; k <3; k++) {
      raw_data = field_reader.read_var_slice("votemper", 0, k);
      for (int i = 0; i<field_view.shape(0); ++i) {
        if (ghost(i)) continue;
        int orca_i = ij(i, 0);
        int orca_j = ij(i, 1);
        size_t raw_idx = glbarray_offset + orca_j * glbarray_jstride + orca_i;
        if (field_view(i, k) != raw_data[raw_idx]) {
            std::cout << "mismatch: "
                      << " ij(" << i << ", 0) " << ij(i, 0)
                      << " ij(" << i << ", 1) " << ij(i, 1)
                      << " checkerboard " << field_view(i, k) << " "
                      << " raw " << raw_data[raw_idx] << std::endl;
        }
        EXPECT_EQUAL(field_view(i, k), raw_data[raw_idx]);
      }
    }
  }
  SECTION("depth field") {
    atlas::Field field(funcSpace.createField<double>(
                          atlas::option::name("depth") |
                          atlas::option::levels(3)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    field_reader.read_vertical_var("nav_lev", mesh, field_view);
    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    std::vector<double> levels{0,10,100};
    for (int k =0; k <3; k++) {
      for (int i = 0; i<field_view.shape(0); ++i) {
        if (ghost(i)) continue;
        if (levels[k] != field_view(i, k)) {
            std::cout << "mismatch: "
                      << " ij(" << i << ", 0) " << ij(i, 0)
                      << " ij(" << i << ", 1) " << ij(i, 1)
                      << " 1 proc " << levels[k] << " "
                      << " with mesh " << field_view(i, k) << std::endl;
        }
        EXPECT_EQUAL(levels[k], field_view(i, k));
      }
    }
  }
}
*/

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
