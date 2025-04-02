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

#include "orca-jedi/nemo_io/WriteServer.h"
#include "orca-jedi/nemo_io/NemoFieldReader.h"
#include "orca-jedi/nemo_io/AtlasIndex.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

/// \brief helper function to apply a lambda across all MPI ranks in turn when a partition is
///        distributed, and just the root partition if not distributed.
/// \param functor           lambda function to apply.
/// \param partitioner_name  name of the partitioner.
template<typename Functor>
void applyMPISerialised(const Functor& functor, const std::string& partitioner_name) {
  if (partitioner_name == "serial") {
    if (atlas::mpi::comm().rank() == 0) {
      functor();
    }
  } else {
    atlas::mpi::comm().barrier();
    size_t rank = 0;
    while (rank < atlas::mpi::comm().size()) {
      if (rank == atlas::mpi::comm().rank()) {
        functor();
      }
      rank++;
      atlas::mpi::comm().barrier();
    }
  }
}

//-----------------------------------------------------------------------------

CASE("test MPI distributed field array view write to disk") {
  auto partitioner_names = std::vector<std::string>{"serial", "checkerboard"};
  auto grid_names = std::vector<std::string>{"ORCA2_T",
     "../Data/amm1_atlas_grid_spec.yaml"};
  size_t nparts = atlas::mpi::size();
  if (nparts == 1)
    partitioner_names = std::vector<std::string>{"serial"};
  for (std::string grid_name : grid_names) {
    for (std::string partitioner_name : partitioner_names) {
      // setup atlas data before write
      eckit::PathName test_data_path;
      eckit::PathName grid_spec_path(grid_name);
      atlas::Grid grid;
      if (grid_spec_path.exists()) {
        grid = atlas::Grid{atlas::Grid::Spec{grid_spec_path}};
        test_data_path = eckit::PathName{std::string("../testoutput/write_" + grid.type() + "_")
          + partitioner_name + "_" + std::to_string(nparts) + ".nc"};
      } else {
        grid = atlas::OrcaGrid{grid_name};
        test_data_path = eckit::PathName{std::string("../testoutput/write_" + grid_name + "_")
            + partitioner_name + "_" + std::to_string(nparts) + ".nc"};
      }

      auto meshgen_config = grid.meshgenerator();
      atlas::MeshGenerator meshgen(meshgen_config);
      auto partitioner_config = grid.partitioner();
      partitioner_config.set("type", partitioner_name);
      auto partitioner = atlas::grid::Partitioner(partitioner_config);
      auto mesh = meshgen.generate(grid, partitioner);

      auto funcSpace = atlas::functionspace::NodeColumns(mesh);
      std::unique_ptr<AtlasIndexToBufferIndex> atlas2buffer(
        AtlasIndexToBufferIndexCreator::create_unique(grid.type(), mesh));

      if (grid_name == "ORCA2_T") {
        SECTION(partitioner_name + "_" + std::to_string(nparts) + " test ORCA indexing") {
          EXPECT(atlas2buffer->nx() == 182);
          EXPECT(atlas2buffer->ny() == 149);
        }
      }

      auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());

      auto field_ice = funcSpace.createField<double>(
        atlas::option::name("iiceconc") | atlas::option::levels(1));
      auto field_t0_temp = funcSpace.createField<float>(
        atlas::option::name("votemper0") | atlas::option::levels(3));
      auto field_t1_temp = funcSpace.createField<float>(
        atlas::option::name("votemper1") | atlas::option::levels(3));

      field_ice.metadata().set("missing_value", 0);
      field_ice.metadata().set("missing_value_type", "approximately-equals");
      field_ice.metadata().set("missing_value_epsilon", 1e-6);
      field_t0_temp.metadata().set("missing_value", 0);
      field_t0_temp.metadata().set("missing_value_type", "approximately-equals");
      field_t0_temp.metadata().set("missing_value_epsilon", 1e-6);
      field_t1_temp.metadata().set("missing_value", 0);
      field_t1_temp.metadata().set("missing_value_type", "approximately-equals");
      field_t1_temp.metadata().set("missing_value_epsilon", 1e-6);

      atlas::field::MissingValue ice_mv(field_ice);
      atlas::field::MissingValue temp0_mv(field_t0_temp);
      atlas::field::MissingValue temp1_mv(field_t1_temp);

      auto field_view_ice = atlas::array::make_view<double, 2>(field_ice);
      auto field_view_t0_temp = atlas::array::make_view<float, 2>(field_t0_temp);
      auto field_view_t1_temp = atlas::array::make_view<float, 2>(field_t1_temp);

      // fill fields
      {
        const size_t num_nodes = field_view_ice.shape(0);
        for (size_t inode = 0; inode < num_nodes; ++inode) {
          if (ghost(inode)) continue;
          field_view_ice(inode, 0) = inode;
          for (size_t ilevel = 0; ilevel < 3; ++ilevel) {
            field_view_t0_temp(inode, ilevel) = inode + ilevel*5;
            field_view_t1_temp(inode, ilevel) = inode + ilevel*10;
          }
        }
      }

      funcSpace.haloExchange(field_ice);
      funcSpace.haloExchange(field_t0_temp);
      funcSpace.haloExchange(field_t1_temp);

      SECTION(partitioner_name + "_" + std::to_string(nparts)
               + " write data to file via write server") {
        std::vector<util::DateTime> datetimes{util::DateTime("1970-01-01T00:00:00Z"),
                                              util::DateTime("1970-01-02T00:00:00Z")};
        std::shared_ptr<eckit::Timer> eckit_timer =
            std::make_shared<eckit::Timer>("write_server tests: ", oops::Log::debug());
        const bool serial_distribution = (atlas::mpi::size() == 1 || partitioner_name == "serial");
        WriteServer writer(eckit_timer, test_data_path, mesh,
                           datetimes, {1, 2, 3}, serial_distribution);

        writer.write_surf_var<double>("iiceconc", 0, ice_mv, field_view_ice);
        writer.write_vol_var<float>("votemper", 0, temp0_mv, field_view_t0_temp);
        writer.write_vol_var<float>("votemper", 1, temp1_mv, field_view_t1_temp);
      }

      SECTION(partitioner_name + "_" + std::to_string(nparts) + " File exists") {
        auto check_file_exists = [&](){
          // wait up to 20 seconds for the file system...
          for (int wait_count=0; wait_count < 10; ++wait_count) {
            if (test_data_path.exists()) break;
            sleep(2);
          }
          EXPECT(test_data_path.exists());
        };
        applyMPISerialised(check_file_exists, partitioner_name);
      }

      SECTION(partitioner_name + "_" + std::to_string(nparts) + " check latitude/longitude data") {
        auto check_lon_lat = [&](){
          NemoFieldReader field_reader(test_data_path);
          std::vector<atlas::PointXY> data = field_reader.read_locs();

          auto lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());
          const size_t num_nodes = mesh.nodes().size();
          EXPECT(num_nodes <= data.size());
          for (size_t inode = 0; inode < num_nodes; ++inode) {
            if (ghost(inode)) continue;
            const int64_t ibuf = (*atlas2buffer)(inode);
            EXPECT(data[ibuf](0) == lonlat(inode, 0));
            EXPECT(data[ibuf](1) == lonlat(inode, 1));
          }
        };
        applyMPISerialised(check_lon_lat, partitioner_name);
      }

      SECTION(partitioner_name + "_" + std::to_string(nparts) + " check surface data") {
        auto check_surface_data = [&](){
          NemoFieldReader field_reader(test_data_path);
          std::vector<double> data = field_reader.read_var_slice<double>("iiceconc", 0, 0);

          const size_t num_nodes = mesh.nodes().size();
          EXPECT(num_nodes <= data.size());
          for (size_t inode = 0; inode < num_nodes; ++inode) {
            if (ghost(inode)) continue;
            const int64_t ibuf = (*atlas2buffer)(inode);
            if (ice_mv(field_view_ice(inode, 0))) {
              EXPECT(std::abs(data[ibuf] - 1e+20) < 1e-6);
            } else {
              EXPECT(data[ibuf] == field_view_ice(inode, 0));
            }
          }
        };
        applyMPISerialised(check_surface_data, partitioner_name);
      }

      SECTION(partitioner_name + "_" + std::to_string(nparts) + " check volume data") {
        auto check_volume_data = [&](){
          NemoFieldReader field_reader(test_data_path);
          std::vector<float> data_t0_l0 = field_reader.read_var_slice<float>("votemper", 0, 0);
          std::vector<float> data_t0_l1 = field_reader.read_var_slice<float>("votemper", 0, 1);
          std::vector<float> data_t0_l2 = field_reader.read_var_slice<float>("votemper", 0, 2);
          std::vector<float> data_t1_l0 = field_reader.read_var_slice<float>("votemper", 1, 0);
          std::vector<float> data_t1_l1 = field_reader.read_var_slice<float>("votemper", 1, 1);
          std::vector<float> data_t1_l2 = field_reader.read_var_slice<float>("votemper", 1, 2);

          auto check_slice = [ghost, mesh](const std::vector<float>& test_data,
              const std::unique_ptr<AtlasIndexToBufferIndex>& o2b,
              const atlas::array::ArrayView<float, 2>& original_fv,
              const atlas::field::MissingValue mv,
              const size_t ilev) {
            const size_t num_nodes = mesh.nodes().size();
            EXPECT(num_nodes <= test_data.size());
            for (size_t inode = 0; inode < num_nodes; ++inode) {
              if (ghost(inode)) continue;
              const int64_t ibuf = (*o2b)(inode);
              if (mv(original_fv(inode, ilev))) continue;
              if (std::abs(test_data[ibuf] - original_fv(inode, ilev)) >= 1e-7)
                std::cout << "test_data[" << ibuf << "] " << test_data[ibuf]
                          << " original_fv(" << inode << ", " << ilev << ") "
                          << original_fv(inode, ilev) << std::endl;
              EXPECT(std::abs(test_data[ibuf] - original_fv(inode, ilev)) < 1e-7);
            }
          };
          check_slice(data_t0_l0, atlas2buffer, field_view_t0_temp, temp0_mv, 0);
          check_slice(data_t0_l1, atlas2buffer, field_view_t0_temp, temp0_mv, 1);
          check_slice(data_t0_l2, atlas2buffer, field_view_t0_temp, temp0_mv, 2);
          check_slice(data_t1_l0, atlas2buffer, field_view_t1_temp, temp1_mv, 0);
          check_slice(data_t1_l1, atlas2buffer, field_view_t1_temp, temp1_mv, 1);
          check_slice(data_t1_l2, atlas2buffer, field_view_t1_temp, temp1_mv, 2);
        };
        applyMPISerialised(check_volume_data, partitioner_name);
      }
    }
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
