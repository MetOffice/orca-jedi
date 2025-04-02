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
  size_t nx = grid.nx() + grid.haloWest() + grid.haloEast();
  size_t ny = grid.ny() + grid.haloSouth() + grid.haloNorth();
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

  std::vector<double> ice_buf, temp_buf;
  for (size_t i = 0; i < ice_fv.shape(0); ++i) {
    ice_fv(i, 0) = ij(i, 0);
    ice_buf.emplace_back(ij(i, 0));
  }
  temp_buf.resize(3*ice_buf.size());
  for (size_t i = 0; i < temp_fv.shape(0); ++i) {
    temp_fv(i, 0) = ij(i, 0);
    temp_fv(i, 1) = ij(i, 1) + 100;
    temp_fv(i, 2) = ij(i, 0) + 200;
    temp_buf[i                   ] = ij(i, 0);
    temp_buf[i +   ice_buf.size()] = ij(i, 1) + 100;
    temp_buf[i + 2*ice_buf.size()] = ij(i, 0) + 200;
  }

  int rank = atlas::mpi::rank();
  SECTION("writes fields") {
    if (rank == 0) {
      std::vector<util::DateTime> datetimes{util::DateTime("1970-01-01T00:00:00Z")};
      {
        NemoFieldWriter field_writer(test_data_path, datetimes, nx, ny, {1, 2, 3});

        atlas::array::ArrayView<double, 2> lonlat{atlas::array::make_view<double, 2>(
                                                    mesh.nodes().lonlat())};
        std::vector<double> lons;
        std::vector<double> lats;
        for (int i_node = 0; i_node < lonlat.shape(0); ++i_node) {
          lons.emplace_back(lonlat(i_node, 0));
          lats.emplace_back(lonlat(i_node, 1));
        }

        EXPECT(lons.size() == nx*ny);
        EXPECT(lons.size() == lonlat.shape(0));
        EXPECT(lats.size() == nx*ny);
        EXPECT(lats.size() == lonlat.shape(0));
        field_writer.write_dimensions(lats, lons);
        std::cout << "sizes: " << ice_buf.size() << " =! " << nx*ny << std::endl;
        EXPECT(ice_buf.size() == nx*ny);
        EXPECT(ice_buf.size() == ice_fv.shape(0));
        field_writer.write_surf_var("iiceconc", ice_buf, 0);
        EXPECT(temp_buf.size() == 3*nx*ny);
        EXPECT(temp_buf.size() == temp_fv.size());
        field_writer.write_vol_var("votemper", temp_buf, 0);
      }
    }
  }

  // wait up to 20 seconds for the file system.
  for (int wait_count=0; wait_count < 10; ++wait_count) {
    if (test_data_path.exists()) break;
    sleep(2);
  }
  EXPECT(test_data_path.exists());

  SECTION("surface field matches with data in memory") {
    if (rank == 0) {
      NemoFieldReader field_reader(test_data_path);
      std::vector<double> data = field_reader.read_var_slice<double>("iiceconc", 0, 0);
      EXPECT(data.size() == nx*ny);
      EXPECT(data.size() == ice_fv.shape(0));
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
      EXPECT(data.size() == nx*ny);
      EXPECT(data.size() == temp_fv.shape(0));
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 0));
      }
      data.clear();
      data = field_reader.read_var_slice<double>("votemper", 0, 1);
      EXPECT(data.size() == nx*ny);
      EXPECT(data.size() == temp_fv.shape(0));
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 1));
      }
      data.clear();
      data = field_reader.read_var_slice<double>("votemper", 0, 2);
      EXPECT(data.size() == nx*ny);
      EXPECT(data.size() == temp_fv.shape(0));
      for (size_t iNode = 0; iNode < data.size(); ++iNode) {
        if (ghost(iNode)) continue;
        EXPECT_EQUAL(data[iNode], temp_fv(iNode, 2));
      }
    }
  }
}
}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
