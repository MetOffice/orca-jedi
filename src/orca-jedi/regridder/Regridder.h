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
/// target grid:
///   name: L90
/// \endcode
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
  atlas::FunctionSpace sourceFunctionSpace_;
  atlas::FunctionSpace targetFunctionSpace_;
  atlas::Interpolation interpolation_;
};

}  // namespace orcamodel
