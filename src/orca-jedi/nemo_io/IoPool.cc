/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/IoPool.h"

#include <algorithm>
#include <string>
#include <utility>

#include "atlas/runtime/Exception.h"  // IWYU pragma: keep
#include "oops/util/Logger.h"

namespace orcamodel {

namespace {
/// \brief Decide whether a given parent rank should be an I/O rank, choosing
///        ranks spread as evenly as possible across the parent communicator.
///
/// Rank r is selected if it is one of the n_io evenly spaced picks
/// floor(k * size / n_io) for k in [0, n_io). This keeps I/O ranks on distinct
/// nodes when ranks are packed node-by-node, improving aggregate filesystem
/// bandwidth. Returns the pool index [0, n_io) if selected, or -1 otherwise.
int pool_index(size_t rank, size_t size, size_t n_io) {
  for (size_t k = 0; k < n_io; ++k) {
    if ((k * size) / n_io == rank) {
      return static_cast<int>(k);
    }
  }
  return -1;
}
}  // namespace

IoPool::IoPool(const eckit::mpi::Comm& parent, size_t n_io_ranks,
               size_t nx, size_t ny, const std::string& name) {
  const size_t size = parent.size();
  const size_t rank = parent.rank();
  parent_size_ = size;

  n_io_ranks_ = std::max<size_t>(1, std::min({n_io_ranks, size, ny == 0 ? size : ny}));

  const int idx = pool_index(rank, size, n_io_ranks_);
  const int color = (idx >= 0) ? 0 : 1;

  // split() is collective over the parent communicator; every rank must call it.
  // It registers the resulting sub-communicator (for this rank's colour) under
  // \p name in eckit's global registry, so both the I/O group and the non-I/O
  // group hold a communicator under this name and both must free it later.
  eckit::mpi::Comm& split_comm = parent.split(color, name);
  comm_name_ = name;

  if (idx >= 0) {
    io_rank_ = idx;
    io_comm_ = &split_comm;
    // split() orders new-communicator ranks by parent rank, so the rank within
    // the I/O communicator matches the evenly spaced pool index.
    ATLAS_ASSERT(static_cast<size_t>(io_comm_->rank()) == static_cast<size_t>(idx),
        "IoPool: unexpected I/O sub-communicator rank ordering");
  }

  decomposition_ = std::make_unique<YBandDecomposition>(nx, ny, n_io_ranks_);

  oops::Log::trace() << "orcamodel::IoPool::IoPool selected " << n_io_ranks_
                     << " I/O ranks from " << size << " (this rank io_rank="
                     << io_rank_ << ")" << std::endl;
}

IoPool::IoPool(IoPool&& other) noexcept
    : n_io_ranks_(other.n_io_ranks_),
      parent_size_(other.parent_size_),
      io_rank_(other.io_rank_),
      io_comm_(other.io_comm_),
      comm_name_(std::move(other.comm_name_)),
      decomposition_(std::move(other.decomposition_)) {
  // The moved-from pool must not free the communicator in its destructor.
  other.comm_name_.clear();
  other.io_comm_ = nullptr;
}

IoPool::~IoPool() {
  // Free the split communicator so the name can be reused on a later pool. Only
  // ranks that registered the name (all parent ranks, unless moved-from) delete
  // it, and each does so collectively within its own colour group.
  if (!comm_name_.empty() && eckit::mpi::hasComm(comm_name_)) {
    eckit::mpi::deleteComm(comm_name_);
  }
}

}  // namespace orcamodel
