/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/testing/Test.h"
#include "atlas/parallel/mpi/mpi.h"

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

#include "orca-jedi/nemo_io/NemoFieldWriter.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"


namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test parallel serially distributed write field array views") {
  eckit::PathName test_data_path("../testoutput/orca1_t_test_output.nc");

  atlas::OrcaGrid grid("ORCA1_T");
  auto meshgen_config = grid.meshgenerator();

  atlas::MeshGenerator meshgen(meshgen_config);
  auto partitioner_config = grid.partitioner();
  partitioner_config.set("type", "serial");
  auto partitioner = atlas::grid::Partitioner(partitioner_config);
  auto mesh = meshgen.generate(grid, partitioner);
  auto funcSpace = atlas::functionspace::NodeColumns(mesh);

  atlas::Field ice_field(funcSpace.createField<double>(
                        atlas::option::name("iiceconc") |
                        atlas::option::levels(1)));
  auto ice_fv = atlas::array::make_view<double, 2>(ice_field);

    atlas::Field temp_field(funcSpace.createField<double>(
                          atlas::option::name("votemper") |
                          atlas::option::levels(3)));
    auto temp_fv = atlas::array::make_view<double, 2>(temp_field);

  auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
  auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());

  for (size_t i = 0; i < ice_fv.size(); ++i) {
    if (ghost(i)) continue;
    ice_fv(i, 0) = ij(i, 0);
    temp_fv(i, 0) = ij(i, 0);
    temp_fv(i, 1) = ij(i, 1) + 100;
    temp_fv(i, 2) = ij(i, 0) + 200;
  }

  int rank = atlas::mpi::rank();
  SECTION("writes fields") {
    if (rank == 0) {
      std::vector<util::DateTime> datetimes(1);
      datetimes[0] = util::DateTime("1970-01-01T00:00:00Z");
      {
        NemoFieldWriter field_writer(test_data_path, mesh, datetimes,
            {1, 2, 3});
        field_writer.write_surf_var("iiceconc", ice_fv, 0);
        field_writer.write_vol_var("votemper", temp_fv, 0);
      }
      // wait up to 20 seconds for the file system...
      for (int wait_count=0; wait_count < 10; ++wait_count) {
        if (test_data_path.exists()) break;
        sleep(2);
      }
      EXPECT(test_data_path.exists());
    }
  }

  SECTION("surface field matches with data in memory") {
    if (rank == 0) {
      NemoFieldReader field_reader(test_data_path);
      std::vector<double> data = field_reader.read_var_slice<double>("iiceconc", 0, 0);
      for (atlas::idx_t iNode = 0; iNode < ice_fv.shape(0); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], ice_fv(iNode, 0));
      }
    }
  }

  SECTION("depth field matches with data in memory") {
    if (rank == 0) {
      NemoFieldReader field_reader(test_data_path);
      auto vert_levels = field_reader.read_vertical_var<double>("z", 3);

      EXPECT_EQUAL(vert_levels[0], 1);
      EXPECT_EQUAL(vert_levels[1], 2);
      EXPECT_EQUAL(vert_levels[2], 3);
    }
  }


  SECTION("volume field matches with data in memory") {
    if (rank == 0) {
      NemoFieldReader field_reader(test_data_path);
      std::vector<double> data = field_reader.read_var_slice<double>("votemper", 0, 0);
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 0));
      }
      data.clear();
      data = field_reader.read_var_slice<double>("votemper", 0, 1);
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 1));
      }
      data.clear();
      data = field_reader.read_var_slice<double>("votemper", 0, 2);
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 2));
      }
    }
  }

//  SECTION("volume field") {
//    atlas::Field field(funcSpace.createField<double>(
//                          atlas::option::name("votemper") |
//                          atlas::option::levels(3)));
//    auto field_view = atlas::array::make_view<double, 2>(field);
//
//    field_reader.read_volume_var("votemper", mesh, 0, field_view);
//    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
//
//    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
//    std::vector<double> raw_data;
//    for (int k =0; k <3; k++) {
//      raw_data = field_reader.read_var_slice("votemper", 0, k);
//      for (int i = 0; i < field_view.shape(0); ++i) {
//        if (ghost(i)) continue;
//        if (raw_data[i] != field_view(i, k)) {
//            std::cout << "mismatch: "
//                      << " ij(" << i << ", 0) " << ij(i, 0)
//                      << " ij(" << i << ", 1) " << ij(i, 1)
//                      << " 1 proc " << raw_data[i] << " "
//                      << " with mesh " << field_view(i, k) << std::endl;
//        }
//        EXPECT_EQUAL(raw_data[i], field_view(i, k));
//      }
//    }
//  }
//
//  SECTION("depth field") {
//    atlas::Field field(funcSpace.createField<double>(
//                          atlas::option::name("depth") |
//                          atlas::option::levels(3)));
//    auto field_view = atlas::array::make_view<double, 2>(field);
//
//    field_reader.read_vertical_var("nav_lev", mesh, field_view);
//    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));
//
//    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
//    std::vector<double> levels{0, 10, 100};
//    for (int k =0; k <3; k++) {
//      for (int i = 0; i < field_view.shape(0); ++i) {
//        if (ghost(i)) continue;
//        if (levels[k] != field_view(i, k)) {
//            std::cout << "mismatch: "
//                      << " ij(" << i << ", 0) " << ij(i, 0)
//                      << " ij(" << i << ", 1) " << ij(i, 1)
//                      << " 1 proc " << levels[k] << " "
//                      << " with mesh " << field_view(i, k) << std::endl;
//        }
//        EXPECT_EQUAL(levels[k], field_view(i, k));
//      }
//    }
//  }
}

// for domain distributed MPI...
//  bool write_nemo_field = false;
//  if (write_nemo_field) {
//    mpi::comm().barrier();
//    size_t rank = 0;
//    while (rank < atlas::mpi::size()) {
//      if (rank == mpi::rank()) {
//        std::cout << "write data from rank " << rank << std::endl;
//        eckit::PathName pathname{"nemo_test_out.nc"};
//        NemoFieldWriter field_writer(pathname, src.mesh());
//        field_writer.write_surf_var(nemo_variable_name, iiceconc);
//      }
//      rank++;
//      mpi::comm().barrier();
//    }
//  }

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
