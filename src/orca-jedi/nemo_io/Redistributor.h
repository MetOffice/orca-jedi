/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <cstddef>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "orca-jedi/nemo_io/IoPool.h"

namespace orcamodel {

/// \brief Bridges the compute (atlas) decomposition and the I/O-slab
///        decomposition with a cached, partitioner-agnostic all-to-all plan.
///
/// The plan is computed once from the global buffer index of each local grid
/// point and reused for every field, level and time step. Data volume moved
/// equals (roughly) a single copy of the field, but the exchange is a balanced
/// all-to-all onto the I/O pool rather than an all-to-one gather onto a single
/// root, and it composes with concurrent collective writes.
///
/// The caller must pass *all* local nodes, including halo/ghost nodes: on the
/// ORCA grid the extra halo buffer cells (east/west wrap and the north fold) are
/// only reachable through ghost nodes, so omitting them would leave those cells
/// unwritten. A global buffer index may consequently be contributed by more than
/// one rank (its owner plus every rank holding it as a ghost). Duplicates are
/// harmless - each copy carries the same value, so a slab cell written more than
/// once is a last-wins no-op, matching the original gather-and-sort on root.
class Redistributor {
 public:
  /// \param comm          Parent communicator spanning all compute ranks (the
  ///                      I/O ranks are a subset of this communicator).
  /// \param pool          The I/O pool defining slab ownership.
  /// \param node_index    Atlas node index of each local point, including ghost
  ///                      nodes (required for full halo-cell coverage).
  /// \param global_index  Flattened global buffer index (g = j*nx + i) of each
  ///                      local point; parallel to node_index.
  Redistributor(const eckit::mpi::Comm& comm, const IoPool& pool,
                const std::vector<size_t>& node_index,
                const std::vector<size_t>& global_index);

  /// \brief Number of horizontal points in this rank's slab (0 on non-I/O ranks).
  size_t slab_size() const { return slab_size_; }

  /// \brief Compute -> I/O. Gather one horizontal level from the per-node
  ///        \p local_values (indexed by atlas node) into this I/O rank's
  ///        contiguous \p slab_out (resized to slab_size()).
  template <class T>
  void to_io(const std::vector<T>& local_values, std::vector<T>& slab_out) const;

  /// \brief I/O -> compute. Scatter this I/O rank's contiguous \p slab_in back
  ///        onto the entries of the per-node \p local_values (which the caller
  ///        must pre-size to the number of atlas nodes). All contributing nodes,
  ///        including ghost nodes, are filled from the file.
  template <class T>
  void from_io(const std::vector<T>& slab_in, std::vector<T>& local_values) const;

 private:
  const eckit::mpi::Comm* comm_;
  std::vector<int> send_counts_;
  std::vector<int> send_displs_;
  std::vector<int> recv_counts_;
  std::vector<int> recv_displs_;
  std::vector<size_t> send_node_order_;   ///< send slot -> local atlas node index
  std::vector<size_t> recv_slab_offset_;  ///< recv slot -> offset within slab
  size_t slab_size_ = 0;
};

}  // namespace orcamodel
