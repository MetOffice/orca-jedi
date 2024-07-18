/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <fstream>

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
#include "orca-jedi/nemo_io/AtlasIndex.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {


//-----------------------------------------------------------------------------
//
CASE("test MPI distributed reads structured grid field array view") {
  eckit::PathName test_data_path("../Data/amm1_nemo.nc");
  eckit::PathName grid_spec_path("../Data/amm1_atlas_grid_spec.yaml");
  auto partitioner_names = std::vector<std::string>{"serial"};
  for (std::string partitioner_name : partitioner_names) {
    atlas::Grid grid{atlas::Grid::Spec{grid_spec_path}};

    auto meshgen_config = grid.meshgenerator();
    std::cout << "meshgen_config: " << meshgen_config << std::endl;
    atlas::MeshGenerator meshgen(meshgen_config);

    auto partitioner_config = grid.partitioner();
    partitioner_config.set("type", partitioner_name);
    auto partitioner = atlas::grid::Partitioner(partitioner_config);

    auto mesh = meshgen.generate(grid, partitioner);
    std::unique_ptr<AtlasIndexToBufferIndex> atlas_index(
      AtlasIndexToBufferIndexCreator::create_unique(grid.type(), mesh));

    auto funcSpace = atlas::functionspace::NodeColumns(mesh);

    auto eckit_timer = std::make_shared<eckit::Timer>(
      "Geometry(ORCA): ", oops::Log::trace());
    eckit_timer->start();
    ReadServer read_server(eckit_timer, test_data_path, mesh);

    constexpr float kgo_missing_value = -32768.;
    SECTION(partitioner_name + " get field _FillValue") {
      auto missing_value = read_server.read_fillvalue<float>("sossheig");
      if (std::abs(missing_value - kgo_missing_value) >= 1e-6) {
        std::cout << "missing_value: " << missing_value << " == " << kgo_missing_value
            << " diff " << std::abs(missing_value - kgo_missing_value) << std::endl;
      }
      EXPECT(std::abs(missing_value - kgo_missing_value) < 1e-6);

      auto default_missing_value = read_server.read_fillvalue<float>("nav_lev");
      if (std::abs(default_missing_value - std::numeric_limits<float>::lowest()) >= 1e-6) {
        std::cout << "default_missing_value: " << missing_value << " == "
            << std::numeric_limits<float>::lowest()
            << " diff " << std::abs(default_missing_value - std::numeric_limits<float>::lowest())
            << std::endl;
      }
      EXPECT(std::abs(default_missing_value - std::numeric_limits<float>::lowest()) < 1e-6);
    }

    SECTION(partitioner_name + " surface field") {
      std::string nemo_name{"sossheig"};
      atlas::Field field(funcSpace.createField<float>(
            atlas::option::name(nemo_name) |
            atlas::option::levels(1)));
      auto field_view = atlas::array::make_view<float, 2>(field);
      field.metadata().set("missing_value", kgo_missing_value);
      field.metadata().set("missing_value_type", "approximately-equals");
      field.metadata().set("missing_value_epsilon", 1e-6);

      read_server.read_var<float>(nemo_name, 0, field_view);

      std::vector<float> raw_data;
      {
        NemoFieldReader field_reader(test_data_path);
        raw_data = field_reader.read_var_slice<float>(nemo_name, 0, 0);
      }

      auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
      for (size_t i = 0; i < field_view.size(); ++i) {
        if (ghost(i)) continue;
        float raw_value = raw_data[(*atlas_index)(i)];
        if (raw_value != field_view(i, 0)) {
          std::cout << "mismatch: i " << i
            << " " << raw_value << " != " << field_view(i, 0) << std::endl;
        }
        EXPECT_EQUAL(raw_value, field_view(i, 0));
      }
    }

    SECTION(partitioner_name + " volume field") {
      atlas::Field field(funcSpace.createField<double>(
            atlas::option::name("votemper") |
            atlas::option::levels(3)));
      auto field_view = atlas::array::make_view<double, 2>(field);

      NemoFieldReader field_reader(test_data_path);
      read_server.read_var<double>("votemper", 0, field_view);

      auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
      std::vector<double> raw_data;
      for (int k =0; k <3; k++) {
        raw_data = field_reader.read_var_slice<double>("votemper", 0, k);
        for (int i = 0; i < field_view.shape(0); ++i) {
          double raw_value = raw_data[(*atlas_index)(i)];
          if (ghost(i)) continue;
          if (raw_value != field_view(i, k)) {
            std::cout << "mismatch: "
              << raw_value << " "
              << " with mesh " << field_view(i, k) << std::endl;
          }
          EXPECT_EQUAL(raw_value, field_view(i, k));
        }
      }
    }
  }
}

CASE("test MPI distributed reads orca grid field array view") {
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
    std::unique_ptr<AtlasIndexToBufferIndex> orca_index(
      AtlasIndexToBufferIndexCreator::create_unique(grid.type(), mesh));

    std::vector<double> raw_data;
    {
      NemoFieldReader field_reader(test_data_path);
      raw_data = field_reader.read_var_slice<double>("iiceconc", 0, 0);
    }

    atlas::Field field(funcSpace.createField<double>(
          atlas::option::name("iiceconc") |
          atlas::option::levels(1)));
    auto field_view = atlas::array::make_view<double, 2>(field);

    auto eckit_timer = std::make_shared<eckit::Timer>(
      "Geometry(ORCA): ", oops::Log::trace());
    eckit_timer->start();
    ReadServer read_server(eckit_timer, test_data_path, mesh);
    read_server.read_var<double>("iiceconc", 0, field_view);

    auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
    SECTION(partitioner_name + " surface field") {
      for (size_t i = 0; i < field_view.size(); ++i) {
        if (ghost(i)) continue;
        double raw_value = raw_data[(*orca_index)(ij(i, 0), ij(i, 1))];
        if (raw_value != field_view(i, 0)) {
          std::cout << "mismatch: "
            << " ij(" << i << ", 0) " << ij(i, 0)
            << " ij(" << i << ", 1) " << ij(i, 1)
            << " 1 proc " << raw_value << " "
            << " with mesh " << field_view(i, 0) << std::endl;
        }
        EXPECT_EQUAL(raw_value, field_view(i, 0));
      }
    }

    SECTION(partitioner_name + " volume field") {
      atlas::Field field(funcSpace.createField<double>(
            atlas::option::name("votemper") |
            atlas::option::levels(3)));
      auto field_view = atlas::array::make_view<double, 2>(field);

      NemoFieldReader field_reader(test_data_path);
      read_server.read_var<double>("votemper", 0, field_view);

      auto ij = atlas::array::make_view<int32_t, 2>(mesh.nodes().field("ij"));

      auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());
      std::vector<double> raw_data;
      for (int k =0; k <3; k++) {
        raw_data = field_reader.read_var_slice<double>("votemper", 0, k);
        for (int i = 0; i < field_view.shape(0); ++i) {
          double raw_value = raw_data[(*orca_index)(ij(i, 0), ij(i, 1))];
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

      read_server.read_vertical_var<double>("nav_lev", field_view);

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
