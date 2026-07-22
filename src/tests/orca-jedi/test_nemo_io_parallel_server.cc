/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "eckit/log/Timer.h"
#include "eckit/filesystem/PathName.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"  // IWYU pragma: keep
#include "atlas/util/Config.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/mesh.h"  // IWYU pragma: keep
#include "atlas/grid.h"  // IWYU pragma: keep
#include "atlas/meshgenerator.h"  // IWYU pragma: keep
#include "atlas/field/Field.h"

#include "atlas-orca/grid/OrcaGrid.h"

#include "orca-jedi/nemo_io/WriteServer.h"
#include "orca-jedi/nemo_io/ReadServer.h"
#include "orca-jedi/nemo_io/AtlasIndex.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("parallel WriteServer + ReadServer round-trip") {
  const size_t nparts = atlas::mpi::size();
  auto partitioner_names = std::vector<std::string>{"serial", "checkerboard"};
  if (nparts == 1) {
    partitioner_names = std::vector<std::string>{"serial"};
  }

  // Values are keyed on the global buffer index so that a node holds the same
  // value regardless of which rank owns it. This makes the round-trip exact.
  const auto surf_val = [](int64_t g) { return static_cast<double>(g) * 0.5 + 1.0; };
  const auto vol_val = [](int64_t g, size_t lev) {
    return static_cast<float>(static_cast<double>(g) * 0.25
                              + static_cast<double>(lev) * 1000.0 + 2.0);
  };

  for (const std::string& partitioner_name : partitioner_names) {
    atlas::OrcaGrid grid("ORCA2_T");
    auto meshgen = atlas::MeshGenerator(grid.meshgenerator());
    auto partitioner_config = grid.partitioner();
    partitioner_config.set("type", partitioner_name);
    auto partitioner = atlas::grid::Partitioner(partitioner_config);
    auto mesh = meshgen.generate(grid, partitioner);
    auto funcSpace = atlas::functionspace::NodeColumns(mesh);

    std::unique_ptr<AtlasIndexToBufferIndex> atlas2buffer(
        AtlasIndexToBufferIndexCreator::create_unique(grid.type(), mesh));
    auto ghost = atlas::array::make_view<int32_t, 1>(mesh.nodes().ghost());

    // Fields without missing_value metadata so values pass through untouched.
    auto field_ice = funcSpace.createField<double>(
        atlas::option::name("iiceconc") | atlas::option::levels(1));
    auto field_temp = funcSpace.createField<float>(
        atlas::option::name("votemper") | atlas::option::levels(3));
    atlas::field::MissingValue ice_mv(field_ice);
    atlas::field::MissingValue temp_mv(field_temp);

    auto view_ice = atlas::array::make_view<double, 2>(field_ice);
    auto view_temp = atlas::array::make_view<float, 2>(field_temp);

    const size_t num_nodes = mesh.nodes().size();
    for (size_t inode = 0; inode < num_nodes; ++inode) {
      if (ghost(inode)) continue;
      const int64_t g = (*atlas2buffer)(inode);
      view_ice(inode, 0) = surf_val(g);
      for (size_t lev = 0; lev < 3; ++lev) {
        view_temp(inode, lev) = vol_val(g, lev);
      }
    }
    funcSpace.haloExchange(field_ice);
    funcSpace.haloExchange(field_temp);

    std::vector<size_t> io_rank_options{0};   // 0 => every rank is an I/O rank
    if (nparts > 1) io_rank_options.push_back(1);

    for (size_t n_io : io_rank_options) {
      SECTION(partitioner_name + "_" + std::to_string(nparts)
              + " n_io=" + std::to_string(n_io)) {
        const eckit::PathName path(std::string("../testoutput/parallel_roundtrip_")
            + partitioner_name + "_" + std::to_string(nparts)
            + "_nio" + std::to_string(n_io) + ".nc");
        const std::vector<util::DateTime> datetimes{
            util::DateTime("1970-01-01T00:00:00Z")};

        // --- parallel write (scoped so the backend closes the file) ---
        {
          std::shared_ptr<eckit::Timer> timer = std::make_shared<eckit::Timer>(
              "parallel_server write: ", oops::Log::debug());
          WriteServer writer(timer, path, mesh, datetimes, {1.0, 2.0, 3.0},
                             /*is_serial=*/false, /*parallel_io=*/true, n_io);
          writer.write_surf_var<double>("iiceconc", 0, ice_mv, view_ice);
          writer.write_vol_var<float>("votemper", 0, temp_mv, view_temp);
        }
        atlas::mpi::comm().barrier();

        // --- parallel read-back into fresh fields ---
        auto rt_ice = funcSpace.createField<double>(
            atlas::option::name("iiceconc") | atlas::option::levels(1));
        auto rt_temp = funcSpace.createField<float>(
            atlas::option::name("votemper") | atlas::option::levels(3));
        auto rt_view_ice = atlas::array::make_view<double, 2>(rt_ice);
        auto rt_view_temp = atlas::array::make_view<float, 2>(rt_temp);

        {
          std::shared_ptr<eckit::Timer> timer = std::make_shared<eckit::Timer>(
              "parallel_server read: ", oops::Log::debug());
          ReadServer reader(timer, path, mesh, /*parallel_io=*/true, n_io);
          reader.read_var<double>("iiceconc", 0, rt_view_ice);
          reader.read_var<float>("votemper", 0, rt_view_temp);
        }

        // Compare on owned (non-ghost) nodes: the round-trip must be exact.
        for (size_t inode = 0; inode < num_nodes; ++inode) {
          if (ghost(inode)) continue;
          EXPECT(rt_view_ice(inode, 0) == view_ice(inode, 0));
          for (size_t lev = 0; lev < 3; ++lev) {
            EXPECT(rt_view_temp(inode, lev) == view_temp(inode, lev));
          }
        }

        atlas::mpi::comm().barrier();
        if (atlas::mpi::rank() == 0) {
          std::remove(path.fullName().asString().c_str());
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
  return orcamodel::test::run(argc, argv);
}
