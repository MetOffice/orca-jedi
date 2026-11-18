/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <memory>
#include <string>

#include "eckit/mpi/Comm.h"

#include "orca-jedi/nemo_io/SlabDecomposition.h"

namespace orcamodel {

/// \brief A subset of MPI ranks dedicated to file I/O, together with the
///        decomposition of the global grid across those ranks.
///
/// The pool selects a configurable number of I/O ranks spread evenly across the
/// parent communicator (to avoid clustering them on a single compute node) and
/// splits off a sub-communicator containing only those ranks. The file is later
/// opened collectively on that sub-communicator, so non-I/O ranks never touch
/// the filesystem. The number of I/O ranks is a runtime knob: it should
/// typically be tuned towards the number of Lustre stripes / OSTs on the target
/// filesystem rather than the total rank count.
class IoPool {
 public:
  /// \param parent        Communicator to draw I/O ranks from (usually the model
  ///                      communicator).
  /// \param n_io_ranks    Desired number of I/O ranks (clamped to [1, parent
  ///                      size] and to the number of grid rows).
  /// \param nx            Global x extent of the grid.
  /// \param ny            Global y extent of the grid.
  /// \param name          Unique name used to register the split communicator.
  IoPool(const eckit::mpi::Comm& parent, size_t n_io_ranks,
         size_t nx, size_t ny, const std::string& name = "orca_io_pool");

  /// \brief Frees the split communicator from the global eckit MPI registry so
  ///        the pool name can be reused (e.g. across successive writes).
  ~IoPool();

  IoPool(IoPool&& other) noexcept;
  IoPool(const IoPool&) = delete;
  IoPool& operator=(IoPool&&) = delete;
  IoPool& operator=(const IoPool&) = delete;

  /// \brief Whether this rank participates in I/O.
  bool is_io_rank() const { return io_rank_ >= 0; }

  /// \brief This rank's index within the I/O pool, or -1 if not an I/O rank.
  int io_rank() const { return io_rank_; }

  /// \brief Number of ranks in the pool.
  size_t n_io_ranks() const { return n_io_ranks_; }

  /// \brief The parent-communicator rank of a given pool member. Mirrors the
  ///        even-spacing selection used in the constructor.
  size_t parent_rank_of(size_t pool_rank) const {
    return (pool_rank * parent_size_) / n_io_ranks_;
  }

  /// \brief The I/O sub-communicator (only valid to use collectively on I/O
  ///        ranks).
  const eckit::mpi::Comm& io_comm() const { return *io_comm_; }

  /// \brief The grid decomposition across the pool.
  const SlabDecomposition& decomposition() const { return *decomposition_; }

 private:
  size_t n_io_ranks_;
  size_t parent_size_ = 1;
  int io_rank_ = -1;
  const eckit::mpi::Comm* io_comm_ = nullptr;
  std::string comm_name_;
  std::unique_ptr<SlabDecomposition> decomposition_;
};

}  // namespace orcamodel
