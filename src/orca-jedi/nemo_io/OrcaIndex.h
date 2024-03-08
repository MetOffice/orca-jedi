/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <algorithm>
#include <sstream>
#include <limits>
#include <map>
#include <utility>
#include <string>
#include <memory>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "atlas/mesh.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

/// \brief Indexer from the orca netcdf i, j field to an index of a 1D buffer
///        of orca data.
struct OrcaIndexToBufferIndex {
 private:
  atlas::OrcaGrid orcaGrid_;
  atlas::Mesh mesh_;
  int32_t ix_glb_max;
  int32_t iy_glb_max;
  int32_t glbarray_offset;
  int32_t glbarray_jstride;
  size_t nx_;
  size_t ny_;

 public:
  size_t nx() const {return nx_;}
  size_t ny() const {return ny_;}

  explicit OrcaIndexToBufferIndex(const atlas::Mesh& mesh) : orcaGrid_(mesh.grid()), mesh_(mesh) {
    iy_glb_max = orcaGrid_.ny() + orcaGrid_.haloNorth() - 1;
    ix_glb_max = orcaGrid_.nx() + orcaGrid_.haloEast() - 1;

    nx_ = orcaGrid_.nx() + orcaGrid_.haloWest() + orcaGrid_.haloEast();
    ny_ = orcaGrid_.ny() + orcaGrid_.haloSouth() + orcaGrid_.haloNorth();

    // vector of local indices: necessary for remote indices of ghost nodes
    int iy_glb_min = -orcaGrid_.haloSouth();
    int ix_glb_min = -orcaGrid_.haloWest();
    glbarray_offset  = -(nx_ * iy_glb_min) - ix_glb_min;
    glbarray_jstride = nx_;
  }

  /// \brief Index of a 1D array corresponding to point i, j
  /// \param i
  /// \param j
  /// \return index of a matching 1D array
  int64_t operator()(const int i, const int j) const {
      ATLAS_ASSERT(i <= ix_glb_max,
          std::to_string(i) + " > " + std::to_string(ix_glb_max));
      ATLAS_ASSERT(j <= iy_glb_max,
          std::to_string(j) + " > " + std::to_string(iy_glb_max));
      return glbarray_offset + j * glbarray_jstride + i;
  }

  /// \brief buffer Index of buffer index corresponding to mesh inode
  /// \param inode
  /// \return index of a matching 1D array
  int64_t operator()(const size_t inode) const {
    auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));
    const size_t i = ij(inode, 0) + orcaGrid_.haloWest();
    const size_t j = ij(inode, 1) + orcaGrid_.haloSouth();
    return j*nx_ + i;
  }
};

}  // namespace orcamodel
