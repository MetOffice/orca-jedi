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
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "atlas/mesh.h"
#include "atlas/grid.h"
#include "atlas/parallel/omp/omp.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

/// \brief Interface from netcdf i, j coordinates in a field to 1D index of atlas
///        field data.
class AtlasIndexToBufferIndex {
 public:
  virtual ~AtlasIndexToBufferIndex() {}
  virtual int64_t operator()(const int i, const int j) const = 0;
  virtual int64_t operator()(const size_t inode) const = 0;
  virtual std::pair<int, int> ij(const size_t inode) const  = 0;
  virtual size_t nx() const = 0;
  virtual size_t ny() const = 0;
};

/// \brief Indexer from the orca netcdf i, j field to an index of a 1D buffer
///        of orca data.
class OrcaIndexToBufferIndex : public AtlasIndexToBufferIndex {
 private:
  atlas::OrcaGrid orcaGrid_;
  int32_t ix_glb_max;
  int32_t iy_glb_max;
  int32_t glbarray_offset;
  int32_t glbarray_jstride;
  size_t nx_;
  size_t ny_;
  const atlas::Mesh mesh_;

 public:
  static std::string name() {return "ORCA";}
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
      ATLAS_ASSERT(i >= ix_glb_min,
          std::to_string(i) + " < " + std::to_string(ix_glb_min));
      ATLAS_ASSERT(j >= iy_glb_min,
          std::to_string(j) + " < " + std::to_string(iy_glb_min));
      return glbarray_offset + j * glbarray_jstride + i;
  }

  /// \brief Index of a 1D array corresponding to a mesh node index
  /// \param inode
  /// \return index of a matching 1D array
  int64_t operator()(const size_t inode) const {
      auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));
      return (*this)(ij(inode, 0), ij(inode, 1));
  }

  /// \brief i, j pair corresponding to a node number. Only use this for diagnostic purposes as
  ///        it is not robust at the orca halo.
  /// \param inode
  /// \return std::pair of the 2D indices corresponding to the node.
  std::pair<int, int> ij(const size_t inode) const {
    auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));
    int i = ij(inode, 0) >= 0 ? ij(inode, 0) : nx_ + ij(inode, 0);
    const int ci = i >= nx_ ? i - nx_ : i;
    int j = ij(inode, 1) + 1;
    const int cj = j >= ny_ ? ny_ - 1 : j;
    return std::pair<int, int>{ci, cj};
  }
};

/// \brief Indexer from the regular lon lat structured grid netcdf i, j field to an index of a
///        1D buffer of regular lon lat field data.
class RegLonLatIndexToBufferIndex : public AtlasIndexToBufferIndex {
 public:
  static std::string name() {return "structured";}
  size_t nx() const {return nx_;}
  size_t ny() const {return ny_;}

  explicit RegLonLatIndexToBufferIndex(const atlas::Mesh& mesh) : grid_(mesh.grid()) {
    ATLAS_ASSERT(grid_.regular(), "RegLonLatIndexToBufferIndex only works with regular grids");
    nx_ = grid_.nx(0);
    ny_ = grid_.ny();

    iy_glb_max = ny_ - 1;
    ix_glb_max = nx_ - 1;

    const size_t num_nodes = mesh.nodes().size();
    inode2ij.resize(num_nodes);
    auto gidx = atlas::array::make_view<atlas::gidx_t, 1>(mesh.nodes().global_index());
    atlas_omp_parallel_for(size_t inode = 0; inode < num_nodes; ++inode) {
      const int global_index = gidx(inode) - 1;
      const int i = global_index % nx_;
      const int j = std::floor(global_index / nx_);
      inode2ij[inode] = std::pair<int, int>(i, j);
    }
    ATLAS_ASSERT(inode2ij.size() <= nx_*ny_,
        std::to_string(inode2ij.size()) + " > " + std::to_string(nx_*ny_));
  }

  /// \brief Index of a 1D array corresponding to point i, j
  /// \param i
  /// \param j
  /// \return index of a matching 1D array
  int64_t operator()(const int i, const int j) const {
    ATLAS_ASSERT(i >= 0, std::to_string(i) + " < 0");
    ATLAS_ASSERT(i <= ix_glb_max,
        std::to_string(i) + " > " + std::to_string(ix_glb_max));
    ATLAS_ASSERT(j >= 0, std::to_string(j) + " < 0");
    ATLAS_ASSERT(j <= iy_glb_max,
        std::to_string(j) + " > " + std::to_string(iy_glb_max));
    return j * nx_ + i;
  }
  /// \brief Index of a 1D array corresponding to a mesh node index
  /// \param inode
  /// \return index of a matching 1D array
  int64_t operator()(const size_t inode) const {
    auto[i, j] = inode2ij[inode];
    return (*this)(i, j);
  }

  /// \brief i, j pair corresponding to a node number.
  /// \param inode
  /// \return std::pair of the 2D indices corresponding to the node.
  std::pair<int, int> ij(const size_t inode) const {
    return inode2ij[inode];
  }

 private:
  atlas::StructuredGrid grid_;
  size_t nx_;
  size_t ny_;
  int32_t ix_glb_max;
  int32_t iy_glb_max;
  std::vector<std::pair<int, int>> inode2ij;
};

/// \brief Factory for creating AtlasIndexToBufferIndex objects.
class AtlasIndexToBufferIndexCreator {
 public:
  static void register_type(std::string name, AtlasIndexToBufferIndexCreator *factory) {
    get_factory()[name] = factory;
  }
  virtual std::unique_ptr<AtlasIndexToBufferIndex> create_unique(const atlas::Mesh& mesh) = 0;
  static std::unique_ptr<AtlasIndexToBufferIndex> create_unique(std::string name,
      const atlas::Mesh& mesh) {
    std::unique_ptr<AtlasIndexToBufferIndex> AtlasIndexToBufferIndex =
      std::move(get_factory()[name]->create_unique(mesh));
    return AtlasIndexToBufferIndex;
  }
  static std::map<std::string, AtlasIndexToBufferIndexCreator*> &get_factory() {
    static std::map<std::string, AtlasIndexToBufferIndexCreator*> creator_map;
    return creator_map;
  }
};

/// \brief Factory for creating OrcaIndexToBufferIndex objects.
class OrcaIndexToBufferIndexCreator : public AtlasIndexToBufferIndexCreator {
 public:
  OrcaIndexToBufferIndexCreator()
    { AtlasIndexToBufferIndexCreator::register_type(OrcaIndexToBufferIndex::name(), this); }
  std::unique_ptr<AtlasIndexToBufferIndex> create_unique(const atlas::Mesh& mesh) {
    std::unique_ptr<AtlasIndexToBufferIndex> o2b_index(new OrcaIndexToBufferIndex(mesh));
    return o2b_index;
  }
};
static OrcaIndexToBufferIndexCreator orcaIndexToBufferIndexCreator;

/// \brief Factory for creating RegLonLatIndexToBufferIndex objects.
class RegLonLatIndexToBufferIndexCreator : public AtlasIndexToBufferIndexCreator {
 public:
  RegLonLatIndexToBufferIndexCreator()
    { AtlasIndexToBufferIndexCreator::register_type(RegLonLatIndexToBufferIndex::name(), this); }
  std::unique_ptr<AtlasIndexToBufferIndex> create_unique(const atlas::Mesh& mesh) {
    std::unique_ptr<AtlasIndexToBufferIndex> r2b_index(new RegLonLatIndexToBufferIndex(mesh));
    return r2b_index;
  }
};
static RegLonLatIndexToBufferIndexCreator regLonLatIndexToBufferIndexCreator;

}  // namespace orcamodel
