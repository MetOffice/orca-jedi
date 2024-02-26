/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "eckit/log/Bytes.h"
#include "eckit/testing/Test.h"

#include "oops/util/DateTime.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"
#include "atlas/util/Config.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/field/Field.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "eckit/exception/Exceptions.h"

#include "orca-jedi/nemo_io/ReadServer.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"
#include "orca-jedi/nemo_io/OrcaIndex.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {


//-----------------------------------------------------------------------------

CASE("test MPI distributed reads field array view") {
  eckit::PathName test_data_path("../Data/orca2_t_nemo.nc");

  auto partitioner_names = std::vector<std::string>{"serial", "checkerboard"};
  for (std::string partitioner_name : partitioner_names) {
    atlas::OrcaGrid grid("ORCA2_T");
    auto meshgen_config = grid.meshgenerator();

    atlas::MeshGenerator meshgen(meshgen_config);
    auto partitioner_config = grid.partitioner();
    partitioner_config.set("type", partitioner_name);
    auto partitioner = atlas::grid::Partitioner(partitioner_config);
    auto mesh = meshgen.generate(grid, partitioner);
    auto funcSpace = atlas::functionspace::NodeColumns(mesh);
    auto orca_index = OrcaIndex(mesh.grid());

    std::vector<double> raw_data;
    {
      NemoFieldReader field_reader(test_data_path);
      raw_data = field_reader.read_var_slice("iiceconc", 0, 0);
    }

    atlas::Field field(funcSpace.createField<double>(
          atlas::option::name("iiceconc") |
          atlas::option::levels(1)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    auto eckit_timer = std::make_shared<eckit::Timer>(
      "Geometry(ORCA): ", oops::Log::trace());
    eckit_timer->start();
    ReadServer read_server(eckit_timer, test_data_path, mesh);
    read_server.read_var("iiceconc", 0, field_view);

    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    for (size_t i = 0; i < field_view.size(); ++i) {
      if (ghost(i)) continue;
      double raw_value = raw_data[orca_index(ij(i, 0), ij(i, 1))];
      if (raw_value != field_view(i, 0)) {
        std::cout << "mismatch: "
          << " ij(" << i << ", 0) " << ij(i, 0)
          << " ij(" << i << ", 1) " << ij(i, 1)
          << " 1 proc " << raw_value << " "
          << " with mesh " << field_view(i, 0) << std::endl;
      }
      EXPECT_EQUAL(raw_value, field_view(i, 0));
    }

    SECTION(partitioner_name + " volume field") {
      atlas::Field field(funcSpace.createField<double>(
            atlas::option::name("votemper") |
            atlas::option::levels(3)));
      auto field_view = atlas::array::make_view<double, 2>(field);

      NemoFieldReader field_reader(test_data_path);
      read_server.read_var("votemper", 0, field_view);

      auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

      auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
      std::vector<double> raw_data;
      for (int k =0; k <3; k++) {
        raw_data = field_reader.read_var_slice("votemper", 0, k);
        for (int i = 0; i < field_view.shape(0); ++i) {
          double raw_value = raw_data[orca_index(ij(i, 0), ij(i, 1))];
          if (ghost(i)) continue;
          if (raw_value != field_view(i, k)) {
            std::cout << "mismatch: "
              << " ij(" << i << ", 0) " << ij(i, 0)
              << " ij(" << i << ", 1) " << ij(i, 1)
              << " 1 proc " << raw_value << " "
              << " with mesh " << field_view(i, k) << std::endl;
          }
          EXPECT_EQUAL(raw_value, field_view(i, k));
        }
      }
    }

    SECTION(partitioner_name + " depth field") {
      atlas::Field field(funcSpace.createField<double>(
            atlas::option::name("depth") |
            atlas::option::levels(3)));
      auto field_view = atlas::array::make_view<double, 2>(field);

      read_server.read_vertical_var("nav_lev", field_view);

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
    SECTION(partitioner_name + " get nearest datetime index") {
      util::DateTime test_datetime{"1970-01-01T00:00:00Z"};
      size_t index = read_server.get_nearest_datetime_index(test_datetime);
      EXPECT_EQUAL(index, 0);
    }

    SECTION(partitioner_name + " get field _FillValue") {
      auto missing_value = read_server.read_fillvalue<double>("iiceconc");
      EXPECT_EQUAL(missing_value, -32768.);

      auto default_missing_value = read_server.read_fillvalue<double>("nav_lat");
      EXPECT_EQUAL(default_missing_value, std::numeric_limits<double>::lowest());
    }
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
