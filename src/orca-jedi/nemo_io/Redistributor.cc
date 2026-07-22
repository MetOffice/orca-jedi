/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/Redistributor.h"

#include <numeric>
#include <vector>

#include "atlas/runtime/Exception.h"  // IWYU pragma: keep
#include "oops/util/Logger.h"

namespace orcamodel {

Redistributor::Redistributor(const eckit::mpi::Comm& comm, const IoPool& pool,
                             const std::vector<size_t>& node_index,
                             const std::vector<size_t>& global_index)
    : comm_(&comm) {
  oops::Log::trace() << "orcamodel::Redistributor::Redistributor" << std::endl;
  ATLAS_ASSERT(node_index.size() == global_index.size(),
      "Redistributor: node_index and global_index must be the same length");

  const size_t size = comm.size();
  const size_t n_owned = global_index.size();
  const SlabDecomposition& decomp = pool.decomposition();

  // 1. Count how many owned points go to each destination parent rank.
  send_counts_.assign(size, 0);
  std::vector<size_t> dest(n_owned);
  for (size_t p = 0; p < n_owned; ++p) {
    const size_t owner_pool_rank = decomp.owner(global_index[p]);
    const size_t dest_rank = pool.parent_rank_of(owner_pool_rank);
    ATLAS_ASSERT(dest_rank < size);
    dest[p] = dest_rank;
    ++send_counts_[dest_rank];
  }

  // 2. Send displacements (prefix sum).
  send_displs_.assign(size, 0);
  for (size_t r = 1; r < size; ++r) {
    send_displs_[r] = send_displs_[r - 1] + send_counts_[r - 1];
  }

  // 3. Bucket owned points by destination, recording the source node and the
  //    global index (which travels alongside so receivers can locate each
  //    value within their slab).
  send_node_order_.resize(n_owned);
  std::vector<size_t> send_global(n_owned);
  {
    std::vector<int> cursor = send_displs_;
    for (size_t p = 0; p < n_owned; ++p) {
      const int slot = cursor[dest[p]]++;
      send_node_order_[slot] = node_index[p];
      send_global[slot] = global_index[p];
    }
  }

  // 4. Exchange counts, then build receive displacements.
  recv_counts_.assign(size, 0);
  comm.allToAll(send_counts_, recv_counts_);
  recv_displs_.assign(size, 0);
  for (size_t r = 1; r < size; ++r) {
    recv_displs_[r] = recv_displs_[r - 1] + recv_counts_[r - 1];
  }
  const size_t n_recv = std::accumulate(recv_counts_.begin(), recv_counts_.end(), size_t{0});

  // 5. Exchange the global indices to learn the received layout.
  std::vector<size_t> recv_global(n_recv);
  comm.allToAllv(send_global.data(), send_counts_.data(), send_displs_.data(),
                 recv_global.data(), recv_counts_.data(), recv_displs_.data());

  // 6. On I/O ranks, translate received global indices into slab offsets.
  //
  // The caller feeds *all* local nodes, including halo/ghost nodes, because on
  // the ORCA grid the extra halo buffer cells (east/west wrap and the north
  // fold) are only reachable through ghost nodes. A given global buffer index
  // can therefore arrive from more than one compute rank (the owner plus every
  // rank holding it as a ghost), so the received count is >= the slab size
  // rather than equal to it. Duplicates are harmless: every copy of a halo
  // point carries the same value, so writing a slab cell more than once is a
  // last-wins no-op (this matches the behaviour of the original gather + sort
  // on the root rank). The requirement we *do* enforce is full coverage: every
  // cell of the slab must be written at least once.
  recv_slab_offset_.resize(n_recv);
  if (pool.is_io_rank()) {
    const Hyperslab slab = decomp.slab(pool.io_rank());
    slab_size_ = slab.size();
    const size_t base = slab.global_index_begin();
    std::vector<char> covered(slab_size_, 0);
    for (size_t r = 0; r < n_recv; ++r) {
      ATLAS_ASSERT(recv_global[r] >= base && recv_global[r] < base + slab_size_,
          "Redistributor: received a global index outside this rank's slab");
      const size_t offset = recv_global[r] - base;
      recv_slab_offset_[r] = offset;
      covered[offset] = 1;
    }
    const size_t n_covered =
        std::accumulate(covered.begin(), covered.end(), size_t{0});
    ATLAS_ASSERT(n_covered == slab_size_,
        "Redistributor: slab not fully covered (covered " + std::to_string(n_covered)
        + " of " + std::to_string(slab_size_) + " cells); ensure all nodes,"
        " including ghost nodes, are passed to the Redistributor");
  } else {
    ATLAS_ASSERT(n_recv == 0, "Redistributor: non-I/O rank received data");
  }
}

template <class T>
void Redistributor::to_io(const std::vector<T>& local_values,
                          std::vector<T>& slab_out) const {
  std::vector<T> send_buf(send_node_order_.size());
  for (size_t s = 0; s < send_node_order_.size(); ++s) {
    send_buf[s] = local_values[send_node_order_[s]];
  }

  std::vector<T> recv_buf(recv_slab_offset_.size());
  comm_->allToAllv(send_buf.data(), send_counts_.data(), send_displs_.data(),
                   recv_buf.data(), recv_counts_.data(), recv_displs_.data());

  slab_out.resize(slab_size_);
  for (size_t r = 0; r < recv_slab_offset_.size(); ++r) {
    slab_out[recv_slab_offset_[r]] = recv_buf[r];
  }
}
template void Redistributor::to_io<double>(const std::vector<double>&,
                                           std::vector<double>&) const;
template void Redistributor::to_io<float>(const std::vector<float>&,
                                          std::vector<float>&) const;

template <class T>
void Redistributor::from_io(const std::vector<T>& slab_in,
                            std::vector<T>& local_values) const {
  // Gather from the slab in the same order the values were received during
  // setup, then run the exchange in reverse (send/recv roles swapped).
  std::vector<T> recv_buf(recv_slab_offset_.size());
  for (size_t r = 0; r < recv_slab_offset_.size(); ++r) {
    recv_buf[r] = slab_in[recv_slab_offset_[r]];
  }

  std::vector<T> send_buf(send_node_order_.size());
  comm_->allToAllv(recv_buf.data(), recv_counts_.data(), recv_displs_.data(),
                   send_buf.data(), send_counts_.data(), send_displs_.data());

  for (size_t s = 0; s < send_node_order_.size(); ++s) {
    local_values[send_node_order_[s]] = send_buf[s];
  }
}
template void Redistributor::from_io<double>(const std::vector<double>&,
                                             std::vector<double>&) const;
template void Redistributor::from_io<float>(const std::vector<float>&,
                                            std::vector<float>&) const;

}  // namespace orcamodel
