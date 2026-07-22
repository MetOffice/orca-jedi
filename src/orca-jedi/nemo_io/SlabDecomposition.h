/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <algorithm>
#include <memory>
#include <string>

#include "atlas/runtime/Exception.h"  // IWYU pragma: keep

namespace orcamodel {

/// \brief Description of the contiguous region of the global (ny, nx) horizontal
///        grid owned by a single I/O rank.
///
/// The decomposition splits the grid into horizontal row-bands (whole rows of
/// x), which map onto a single contiguous netCDF hyperslab and a single
/// contiguous run of the flattened global buffer index g = j * nx + i. Levels
/// are handled by the caller looping over level ranges, so this struct is
/// purely two-dimensional and independent of the number of vertical levels.
struct Hyperslab {
  size_t j_begin = 0;   ///< First global row (y index) owned by the I/O rank.
  size_t j_count = 0;   ///< Number of rows owned.
  size_t nx = 0;        ///< Full x extent of the grid (rows are never split).

  /// \brief First flattened global buffer index owned by this slab.
  size_t global_index_begin() const { return j_begin * nx; }

  /// \brief Number of horizontal points owned by this slab.
  size_t size() const { return j_count * nx; }

  /// \brief Whether a flattened global buffer index falls inside this slab.
  bool contains(size_t global_index) const {
    return global_index >= global_index_begin()
        && global_index < global_index_begin() + size();
  }
};

/// \brief Strategy describing how the global horizontal grid is partitioned
///        across the I/O ranks.
///
/// This is deliberately independent of the compute (atlas) decomposition: it
/// only ever reasons about the global flattened buffer index of a grid point.
/// Any MPI partitioner can therefore be paired with any SlabDecomposition; the
/// Redistributor bridges the two. Implementations must guarantee that the slabs
/// tile the whole grid exactly once (no gaps, no overlaps).
class SlabDecomposition {
 public:
  virtual ~SlabDecomposition() = default;

  /// \brief Number of I/O ranks the grid is split across.
  virtual size_t n_io_ranks() const = 0;

  /// \brief The row-band owned by a given I/O rank.
  virtual Hyperslab slab(size_t io_rank) const = 0;

  /// \brief The I/O rank that owns a given flattened global buffer index.
  virtual size_t owner(size_t global_index) const = 0;
};

/// \brief Row-band (split on y) decomposition of the global (ny, nx) grid.
///
/// Rows are distributed as evenly as possible: the first (ny % n_io_ranks)
/// ranks receive one extra row. Each band is fully contiguous both in the file
/// and in the flattened global buffer, which is the friendliest layout for
/// parallel HDF5 collective writes and Lustre striping.
class YBandDecomposition : public SlabDecomposition {
 public:
  YBandDecomposition(size_t nx, size_t ny, size_t n_io_ranks)
      : nx_(nx), ny_(ny),
        n_io_ranks_(std::min(n_io_ranks, ny == 0 ? size_t{1} : ny)),
        base_(ny_ / n_io_ranks_), remainder_(ny_ % n_io_ranks_) {
    ATLAS_ASSERT(n_io_ranks_ >= 1, "YBandDecomposition needs at least one I/O rank");
    ATLAS_ASSERT(nx_ >= 1, "YBandDecomposition needs a non-zero x extent");
  }

  static std::string name() { return "y-band"; }

  size_t n_io_ranks() const override { return n_io_ranks_; }

  Hyperslab slab(size_t io_rank) const override {
    ATLAS_ASSERT(io_rank < n_io_ranks_);
    // Ranks below remainder_ carry one extra row.
    const size_t extra = std::min(io_rank, remainder_);
    const size_t j_begin = io_rank * base_ + extra;
    const size_t j_count = base_ + (io_rank < remainder_ ? 1 : 0);
    return Hyperslab{j_begin, j_count, nx_};
  }

  size_t owner(size_t global_index) const override {
    const size_t j = global_index / nx_;
    // Rows [0, remainder_*(base_+1)) belong to the "fat" ranks.
    const size_t fat_rows = remainder_ * (base_ + 1);
    if (j < fat_rows) {
      return j / (base_ + 1);
    }
    return remainder_ + (j - fat_rows) / base_;
  }

 private:
  size_t nx_;
  size_t ny_;
  size_t n_io_ranks_;
  size_t base_;       ///< Rows per rank before distributing the remainder.
  size_t remainder_;  ///< Number of ranks receiving one extra row.
};

}  // namespace orcamodel
