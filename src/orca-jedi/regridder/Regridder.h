/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/interpolation.h"
#include "atlas/meshgenerator.h"
#include "atlas/mesh.h"
#include "eckit/config/Configuration.h"

namespace orcamodel {

/// \brief Atlas-only regridder class.
///
/// Performs interpolation/regridding between any two atlas FunctionSpaces.
/// No OOPS dependency — designed to be extractable to a standalone tool.
///
/// Configuration example:
/// \code{.yaml}
/// interpolation method:
///   type: finite-element
///   non_linear: missing-if-all-missing
///   cache path: /scratch/regrid/orca025_to_L90   # optional; see below
/// target grid:
///   name: L90
/// \endcode
///
/// \note If "cache path" is set, the interpolation matrix is written to disk on
///   first construction and reloaded on subsequent runs, skipping the (often
///   dominant) matrix-build cost. The stored file is keyed by MPI size and rank
///   because the matrix maps the local source partition to the local target
///   partition, so a cache is only reused for a matching decomposition. Only
///   matrix-based interpolation methods (e.g. finite-element,
///   unstructured-bilinear-lonlat) can be cached.
class Regridder {
 public:
  /// \brief Construct a Regridder from source/target FunctionSpaces and config.
  /// \param conf Configuration specifying interpolation method and options.
  /// \param sourceFunctionSpace The source function space (already constructed).
  /// \param targetFunctionSpace The target function space (already constructed).
  Regridder(const eckit::Configuration& conf,
            const atlas::FunctionSpace& sourceFunctionSpace,
            const atlas::FunctionSpace& targetFunctionSpace);

  /// \brief Execute regridding on a FieldSet.
  /// \param source Source fields to regrid.
  /// \return A new FieldSet on the target function space.
  atlas::FieldSet execute(const atlas::FieldSet& source) const;

  /// \brief Execute regridding on a single Field.
  /// \param source Source field to regrid.
  /// \return A new Field on the target function space.
  atlas::Field execute(const atlas::Field& source) const;

  /// \brief Access the underlying atlas Interpolation object.
  const atlas::Interpolation& interpolation() const { return interpolation_; }

 private:
  /// \brief Build the interpolation object for the given function spaces.
  ///
  /// When the configuration contains a "cache path", the interpolation matrix
  /// is loaded from disk if a matching (per MPI size/rank) cache file exists,
  /// otherwise it is built and written to disk for reuse on later runs.
  static atlas::Interpolation makeInterpolation(
      const eckit::Configuration& conf,
      const atlas::FunctionSpace& source,
      const atlas::FunctionSpace& target);

  atlas::FunctionSpace sourceFunctionSpace_;
  atlas::FunctionSpace targetFunctionSpace_;
  atlas::Interpolation interpolation_;
};

}  // namespace orcamodel
