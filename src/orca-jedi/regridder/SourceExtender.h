/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/mesh.h"

namespace orcamodel {

/// \brief Connectivity type for computing node-to-node adjacency.
enum class AdjacencyType {
  /// Cell-based: neighbours are all nodes sharing a cell with a given node.
  /// On a quad mesh this gives the 8 surrounding nodes (4 edge + 4 corner),
  /// matching the stencil used by NEMOVAR's sim_ext.
  CellBased,
  /// Edge-based: neighbours are nodes connected by a mesh edge.
  /// Typically 4–6 neighbours on a quad mesh. More geometrically principled
  /// but a sparser stencil than NEMOVAR.
  EdgeBased
};

/// \brief Extends source field data into missing-value (land) regions.
///
/// Before regridding from a coarser to a finer grid, target ocean points
/// whose entire interpolation stencil falls on source land will receive
/// missing values. This class fills ("floods") missing-value cells on the
/// source grid by iteratively averaging from valid neighbours, then
/// smooths the newly-flooded cells to remove discontinuities.
///
/// This is analogous to NEMOVAR's sim_ext module (VAR_SRC/SIM/sim_ext.F90).
///
/// Configuration example:
/// \code{.yaml}
/// source extension:
///   flood iterations: 2
///   smooth iterations: 50
///   smooth weight self: 0.35
///   adjacency: cell-based   # or "edge-based"
/// \endcode
class SourceExtender {
 public:
  /// \brief Construct a SourceExtender.
  /// \param mesh The atlas Mesh for the source grid. Connectivity structures
  ///             may be built on the mesh (non-const copy is taken internally
  ///             via atlas's handle semantics).
  /// \param funcSpace The NodeColumns function space (for halo exchange).
  /// \param nFloodIterations Number of flood-fill iterations (cells deep).
  /// \param nSmoothIterations Smoothing iterations on flooded cells.
  /// \param smoothWeightSelf Self-weight in smoothing (NEMOVAR uses 0.35).
  /// \param adjacency Connectivity type for neighbour lookup.
  SourceExtender(const atlas::Mesh& mesh,
                 const atlas::FunctionSpace& funcSpace,
                 int nFloodIterations = 2,
                 int nSmoothIterations = 50,
                 double smoothWeightSelf = 0.35,
                 AdjacencyType adjacency = AdjacencyType::CellBased);

  /// \brief Extend all fields in a FieldSet in-place.
  ///
  /// For each field, missing-value cells adjacent to valid cells are
  /// iteratively filled, then the newly-filled cells are smoothed.
  /// Fields must have "missing_value" metadata set.
  void extend(atlas::FieldSet& fields) const;

  /// \brief Extend a single Field in-place.
  void extend(atlas::Field& field) const;

 private:
  /// Pre-computed node-to-node adjacency list.
  std::vector<std::vector<atlas::idx_t>> adjacency_;

  atlas::FunctionSpace funcSpace_;
  int nFloodIterations_;
  int nSmoothIterations_;
  double smoothWeightSelf_;
};

}  // namespace orcamodel
