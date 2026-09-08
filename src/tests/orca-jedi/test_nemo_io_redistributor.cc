/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "atlas/parallel/mpi/mpi.h"
#include "atlas/array.h"  // IWYU pragma: keep
#include "atlas/grid.h"  // IWYU pragma: keep
#include "atlas/mesh.h"  // IWYU pragma: keep
#include "atlas/meshgenerator.h"  // IWYU pragma: keep

#include "atlas-orca/grid/OrcaGrid.h"

#include "orca-jedi/nemo_io/AtlasIndex.h"
#include "orca-jedi/nemo_io/IoPool.h"
#include "orca-jedi/nemo_io/Redistributor.h"
#include "orca-jedi/nemo_io/SlabDecomposition.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("YBandDecomposition tiles the grid exactly once") {
  const size_t nx = 10;
  const size_t ny = 7;
  for (size_t n_io : std::vector<size_t>{1, 2, 3, 7}) {
    YBandDecomposition decomp(nx, ny, n_io);
    EXPECT(decomp.n_io_ranks() == std::min(n_io, ny));

    // Every global index is owned by exactly the slab that contains it.
    std::vector<int> covered(nx * ny, 0);
    for (size_t r = 0; r < decomp.n_io_ranks(); ++r) {
      const Hyperslab slab = decomp.slab(r);
      for (size_t g = slab.global_index_begin();
           g < slab.global_index_begin() + slab.size(); ++g) {
        EXPECT(decomp.owner(g) == r);
        covered[g] += 1;
      }
    }
    for (size_t g = 0; g < nx * ny; ++g) {
      EXPECT(covered[g] == 1);
    }
  }
}

//-----------------------------------------------------------------------------

CASE("Redistributor compute<->IO round trip") {
  auto partitioner_names = std::vector<std::string>{"serial", "checkerboard"};
  const size_t nparts = atlas::mpi::size();
  if (nparts == 1) {
    partitioner_names = std::vector<std::string>{"serial"};
  }

  const std::string grid_name = "ORCA2_T";
  const atlas::OrcaGrid grid{grid_name};

  for (const std::string& partitioner_name : partitioner_names) {
    auto meshgen_config = grid.meshgenerator();
    atlas::MeshGenerator meshgen(meshgen_config);
    auto partitioner_config = grid.partitioner();
    partitioner_config.set("type", partitioner_name);
    auto partitioner = atlas::grid::Partitioner(partitioner_config);
    auto mesh = meshgen.generate(grid, partitioner);

    std::unique_ptr<AtlasIndexToBufferIndex> atlas2buffer(
        AtlasIndexToBufferIndexCreator::create_unique(grid.type(), mesh));

    const size_t num_nodes = mesh.nodes().size();

    // All local points (including ghost nodes) and their global buffer indices.
    // Ghost nodes are required so that the ORCA halo buffer cells (east/west
    // wrap and the north fold), which no non-ghost node maps to, are covered.
    std::vector<size_t> node_index;
    std::vector<size_t> global_index;
    for (size_t inode = 0; inode < num_nodes; ++inode) {
      node_index.emplace_back(inode);
      global_index.emplace_back(static_cast<size_t>((*atlas2buffer)(inode)));
    }

    // A value that is a known function of the global index, so we can verify
    // placement on the I/O side independently of the compute partition.
    const auto expected = [](size_t g) { return static_cast<double>(g) * 0.5 + 1.0; };

    std::vector<size_t> io_rank_counts{1};
    if (nparts > 1) io_rank_counts.push_back(nparts);

    for (size_t n_io : io_rank_counts) {
      const std::string pool_name =
          "test_pool_" + partitioner_name + "_" + std::to_string(n_io);
      IoPool pool(atlas::mpi::comm(), n_io, atlas2buffer->nx(), atlas2buffer->ny(),
                  pool_name);
      Redistributor redist(atlas::mpi::comm(), pool, node_index, global_index);

      // Fill a per-node buffer for the owned points.
      std::vector<double> local(num_nodes, -999.0);
      for (size_t k = 0; k < node_index.size(); ++k) {
        local[node_index[k]] = expected(global_index[k]);
      }

      // compute -> IO
      std::vector<double> slab;
      redist.to_io(local, slab);

      SECTION(partitioner_name + " n_io=" + std::to_string(n_io) + " slab placement") {
        if (pool.is_io_rank()) {
          const Hyperslab h = pool.decomposition().slab(pool.io_rank());
          EXPECT(slab.size() == h.size());
          for (size_t kk = 0; kk < slab.size(); ++kk) {
            EXPECT(slab[kk] == expected(h.global_index_begin() + kk));
          }
        } else {
          EXPECT(slab.empty());
        }
      }

      // IO -> compute
      SECTION(partitioner_name + " n_io=" + std::to_string(n_io) + " round trip") {
        std::vector<double> local2(num_nodes, -111.0);
        redist.from_io(slab, local2);
        for (size_t k = 0; k < node_index.size(); ++k) {
          EXPECT(local2[node_index[k]] == expected(global_index[k]));
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
