/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/regridder/SourceExtender.h"

#include <set>
#include <vector>

#include "atlas/array/MakeView.h"  // IWYU pragma: keep
#include "atlas/array/DataType.h"  // IWYU pragma: keep
#include "atlas/field/MissingValue.h"
#include "atlas/mesh/Connectivity.h"
#include "atlas/mesh/Nodes.h"
#include "atlas/mesh/actions/BuildEdges.h"
#include "atlas/mesh/actions/BuildNode2CellConnectivity.h"
#include "atlas/parallel/omp/omp.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

namespace orcamodel {

namespace {

/// Build node-to-node adjacency via shared cells (includes diagonal neighbours).
std::vector<std::vector<atlas::idx_t>> buildCellBasedAdjacency(
    atlas::Mesh& mesh) {
  // Ensure node-to-cell connectivity is built
  atlas::mesh::actions::build_node_to_cell_connectivity(mesh);

  const auto& node2cell = mesh.nodes().cell_connectivity();
  const auto& cell2node = mesh.cells().node_connectivity();
  const atlas::idx_t nNodes = mesh.nodes().size();

  std::vector<std::vector<atlas::idx_t>> adj(nNodes);
  // Each iteration builds an independent, pre-allocated adj[jnode] from a
  // thread-local set, reading only const connectivity: safe to parallelise.
  atlas_omp_parallel_for(atlas::idx_t jnode = 0; jnode < nNodes; ++jnode) {
    std::set<atlas::idx_t> neighbours;
    for (atlas::idx_t jc = 0; jc < node2cell.cols(jnode); ++jc) {
      atlas::idx_t icell = node2cell(jnode, jc);
      if (icell < 0) continue;  // missing connectivity entry
      for (atlas::idx_t jn = 0; jn < cell2node.cols(icell); ++jn) {
        atlas::idx_t candidate = cell2node(icell, jn);
        if (candidate != jnode && candidate >= 0) {
          neighbours.insert(candidate);
        }
      }
    }
    adj[jnode].assign(neighbours.begin(), neighbours.end());
  }
  return adj;
}

/// Build node-to-node adjacency via edges (edge-connected neighbours only).
std::vector<std::vector<atlas::idx_t>> buildEdgeBasedAdjacency(
    atlas::Mesh& mesh) {
  // Ensure edges and node-to-edge connectivity are built
  atlas::mesh::actions::build_edges(mesh);
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);

  const auto& node2edge = mesh.nodes().edge_connectivity();
  const auto& edge2node = mesh.edges().node_connectivity();
  const atlas::idx_t nNodes = mesh.nodes().size();

  std::vector<std::vector<atlas::idx_t>> adj(nNodes);
  // Each iteration builds an independent, pre-allocated adj[jnode] from a
  // thread-local set, reading only const connectivity: safe to parallelise.
  atlas_omp_parallel_for(atlas::idx_t jnode = 0; jnode < nNodes; ++jnode) {
    std::set<atlas::idx_t> neighbours;
    for (atlas::idx_t je = 0; je < node2edge.cols(jnode); ++je) {
      atlas::idx_t iedge = node2edge(jnode, je);
      if (iedge < 0) continue;
      atlas::idx_t n_a = edge2node(iedge, 0);
      atlas::idx_t n_b = edge2node(iedge, 1);
      atlas::idx_t neighbour = (n_a == jnode) ? n_b : n_a;
      if (neighbour >= 0) {
        neighbours.insert(neighbour);
      }
    }
    adj[jnode].assign(neighbours.begin(), neighbours.end());
  }
  return adj;
}

}  // namespace

SourceExtender::SourceExtender(const atlas::Mesh& mesh,
                               const atlas::FunctionSpace& funcSpace,
                               int nFloodIterations,
                               int nSmoothIterations,
                               double smoothWeightSelf,
                               AdjacencyType adjacency)
    : funcSpace_(funcSpace),
      nFloodIterations_(nFloodIterations),
      nSmoothIterations_(nSmoothIterations),
      smoothWeightSelf_(smoothWeightSelf) {
  oops::Log::trace() << "orcamodel::SourceExtender constructor starting" << std::endl;
  ASSERT(nFloodIterations >= 0);
  ASSERT(nSmoothIterations >= 0);
  ASSERT(smoothWeightSelf >= 0.0 && smoothWeightSelf <= 1.0);

  // Take a non-const copy via atlas handle semantics (cheap, shared data)
  atlas::Mesh mutableMesh(mesh);

  if (adjacency == AdjacencyType::CellBased) {
    oops::Log::info() << "SourceExtender: building cell-based adjacency"
                      << std::endl;
    adjacency_ = buildCellBasedAdjacency(mutableMesh);
  } else {
    oops::Log::info() << "SourceExtender: building edge-based adjacency"
                      << std::endl;
    adjacency_ = buildEdgeBasedAdjacency(mutableMesh);
  }
  oops::Log::info() << "SourceExtender: adjacency built for "
                    << adjacency_.size() << " nodes" << std::endl;
  oops::Log::trace() << "orcamodel::SourceExtender constructor finished" << std::endl;
}

void SourceExtender::extend(atlas::FieldSet& fields) const {
  oops::Log::trace() << "orcamodel::SourceExtender extend fieldset starting" << std::endl;
  for (atlas::idx_t i = 0; i < fields.size(); ++i) {
    extend(fields[i]);
  }
  oops::Log::trace() << "orcamodel::SourceExtender extend fieldset finished" << std::endl;
}

void SourceExtender::extend(atlas::Field& field) const {
  oops::Log::trace() << "orcamodel::SourceExtender extend field '"
      << field.name() << "' starting" << std::endl;
  if (nFloodIterations_ == 0) {
    oops::Log::trace() << "orcamodel::SourceExtender extend field '"
        << field.name() << "' finished (no flooding)" << std::endl;
    return;
  }

  // Determine missing value from field metadata
  atlas::field::MissingValue mv(field);
  if (!mv) {
    oops::Log::warning()
        << "SourceExtender::extend: field '" << field.name()
        << "' has no missing_value metadata, skipping" << std::endl;
    oops::Log::trace() << "orcamodel::SourceExtender extend field '"
        << field.name() << "' finished (no missing values)" << std::endl;
    return;
  }

  // Dispatch based on field datatype
  if (field.datatype() == atlas::array::DataType::real64()) {
    extendTyped<double>(field);
  } else if (field.datatype() == atlas::array::DataType::real32()) {
    extendTyped<float>(field);
  } else {
    oops::Log::warning()
        << "SourceExtender::extend: field '" << field.name()
        << "' has unsupported datatype (kind="
        << field.datatype().kind() << "), skipping" << std::endl;
    oops::Log::trace() << "orcamodel::SourceExtender extend field '"
        << field.name() << "' finished (unsupported datatype)" << std::endl;
    return;
  }
  oops::Log::trace() << "orcamodel::SourceExtender extend field '"
      << field.name() << "' finished" << std::endl;
}

template <typename T>
void SourceExtender::extendTyped(atlas::Field& field) const {
  atlas::field::MissingValue mv(field);

  const atlas::idx_t nNodes = field.shape(0);
  const atlas::idx_t nLevels = field.levels() > 0 ? field.levels() : 1;

  auto view = atlas::array::make_view<T, 2>(field);
  auto ghost = atlas::array::make_view<int32_t, 1>(
      funcSpace_.ghost());

  // Track which nodes were originally missing (per level).
  // These are the nodes eligible for flooding and smoothing.
  std::vector<std::vector<bool>> originallyMissing(
      nLevels, std::vector<bool>(nNodes, false));
  for (atlas::idx_t k = 0; k < nLevels; ++k) {
    for (atlas::idx_t j = 0; j < nNodes; ++j) {
      originallyMissing[k][j] = mv(view(j, k));
    }
  }

  // Working mask: true = currently valid (ocean). Updated as we flood.
  std::vector<std::vector<bool>> isValid(
      nLevels, std::vector<bool>(nNodes, false));
  for (atlas::idx_t k = 0; k < nLevels; ++k) {
    for (atlas::idx_t j = 0; j < nNodes; ++j) {
      isValid[k][j] = !originallyMissing[k][j];
    }
  }

  // Track which nodes were flooded (for selective smoothing)
  std::vector<std::vector<bool>> wasFlooded(
      nLevels, std::vector<bool>(nNodes, false));

  // --- Flood-fill iterations ---
  oops::Log::debug() << "SourceExtender::extend: flooding field '"
                     << field.name() << "' with " << nFloodIterations_
                     << " flood iterations" << std::endl;
  for (int iter = 0; iter < nFloodIterations_; ++iter) {
    // Snapshot the current validity for this iteration
    auto validSnapshot = isValid;

    for (atlas::idx_t k = 0; k < nLevels; ++k) {
      for (atlas::idx_t j = 0; j < nNodes; ++j) {
        if (ghost(j)) continue;  // ghost nodes filled via haloExchange
        if (validSnapshot[k][j]) continue;  // already valid

        // Check neighbours for valid values
        double weightedSum = 0.0;
        double totalWeight = 0.0;
        for (atlas::idx_t neighbour : adjacency_[j]) {
          if (neighbour >= nNodes) continue;
          if (validSnapshot[k][neighbour]) {
            weightedSum += view(neighbour, k);
            totalWeight += 1.0;
          }
        }

        if (totalWeight > 0.0) {
          // At least one valid neighbour: fill this cell
          view(j, k) = static_cast<T>(weightedSum / totalWeight);
          isValid[k][j] = true;
          wasFlooded[k][j] = true;
        }
      }
    }

    // Halo exchange after each flood iteration
    funcSpace_.haloExchange(field);

    // Re-sync isValid for ghost nodes from the actual field values,
    // since haloExchange may have updated ghost values.
    for (atlas::idx_t k = 0; k < nLevels; ++k) {
      for (atlas::idx_t j = 0; j < nNodes; ++j) {
        if (!ghost(j)) continue;
        bool nowValid = !mv(view(j, k));
        isValid[k][j] = nowValid;
        if (nowValid && originallyMissing[k][j]) {
          wasFlooded[k][j] = true;
        }
      }
    }
  }

  // --- Smoothing on flooded cells only ---
  if (nSmoothIterations_ > 0) {
    const double neighbourWeight = 1.0 - smoothWeightSelf_;
    std::vector<T> tmp(nNodes);
    oops::Log::debug() << "SourceExtender::extend: smoothing field '"
                       << field.name() << "' with " << nSmoothIterations_
                       << " smooth iterations" << std::endl;

    for (int iter = 0; iter < nSmoothIterations_; ++iter) {
      for (atlas::idx_t k = 0; k < nLevels; ++k) {
        // Compute pass writes only tmp[j] (distinct per j); all reads of the
        // field view and the bool masks are read-only, so this is safe to
        // parallelise over nodes.
        atlas_omp_parallel_for(atlas::idx_t j = 0; j < nNodes; ++j) {
          if (ghost(j)) continue;             // ghost nodes via haloExchange
          if (!wasFlooded[k][j]) continue;  // only smooth flooded cells

          double neighbourSum = 0.0;
          int neighbourCount = 0;
          for (atlas::idx_t neighbour : adjacency_[j]) {
            if (neighbour >= nNodes) continue;
            if (!isValid[k][neighbour]) continue;  // skip still-missing cells
            neighbourSum += view(neighbour, k);
            ++neighbourCount;
          }

          if (neighbourCount > 0) {
            tmp[j] = static_cast<T>(smoothWeightSelf_ * view(j, k)
                   + neighbourWeight * neighbourSum / neighbourCount);
          } else {
            tmp[j] = view(j, k);
          }
        }
        // Apply smoothed values: each j writes a distinct view(j, k).
        atlas_omp_parallel_for(atlas::idx_t j = 0; j < nNodes; ++j) {
          if (wasFlooded[k][j]) {
            view(j, k) = tmp[j];
          }
        }
      }

      // Halo exchange after each smoothing iteration
      funcSpace_.haloExchange(field);
    }
  }

  oops::Log::debug() << "SourceExtender::extend: extended field '"
                     << field.name() << "' with " << nFloodIterations_
                     << " flood + " << nSmoothIterations_
                     << " smooth iterations" << std::endl;
}

}  // namespace orcamodel
