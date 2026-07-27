/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/regridder/Regridder.h"

#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/interpolation/Cache.h"
#include "atlas/linalg/sparse/MakeEckitSparseMatrix.h"
#include "atlas/linalg/sparse/MakeSparseMatrixStorageEckit.h"
#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"
#include "eckit/linalg/SparseMatrix.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

namespace orcamodel {

namespace {

/// \brief Per-rank cache filename for the interpolation matrix.
///
/// The matrix maps the local source partition to the local target partition, so
/// it is only valid for the same MPI decomposition (rank count + partitioner)
/// that produced it. Encoding size and rank in the filename prevents an
/// incompatible cache from being reused.
std::string cacheFileName(const std::string& base) {
  const eckit::mpi::Comm& comm = eckit::mpi::comm();
  return base + ".mpi" + std::to_string(comm.size())
              + "_" + std::to_string(comm.rank()) + ".eckit";
}

}  // namespace

atlas::Interpolation Regridder::makeInterpolation(
    const eckit::Configuration& conf,
    const atlas::FunctionSpace& source,
    const atlas::FunctionSpace& target) {
  std::string cachePath;
  conf.get("cache path", cachePath);

  if (cachePath.empty()) {
    return atlas::Interpolation(conf, source, target);
  }

  const eckit::PathName file(cacheFileName(cachePath));

  if (file.exists()) {
    oops::Log::info() << "orcamodel::Regridder loading interpolation matrix "
        << "cache from " << file << std::endl;
    eckit::linalg::SparseMatrix eckitMatrix;
    eckitMatrix.load(file);
    atlas::interpolation::MatrixCache cache(
        atlas::linalg::make_sparse_matrix_storage(std::move(eckitMatrix)));
    return atlas::Interpolation(conf, source, target, cache);
  }

  oops::Log::info() << "orcamodel::Regridder building interpolation matrix; "
      << "will cache to " << file << std::endl;
  atlas::Interpolation interpolation(conf, source, target);

  atlas::interpolation::MatrixCache cache(interpolation.createCache());
  if (cache) {
    eckit::linalg::SparseMatrix eckitMatrix =
        atlas::linalg::make_eckit_sparse_matrix(cache.matrix());
    eckitMatrix.save(file);
    oops::Log::info() << "orcamodel::Regridder cached interpolation matrix to "
        << file << std::endl;
  } else {
    oops::Log::warning() << "orcamodel::Regridder: interpolation method '"
        << conf.getString("type", "<unknown>") << "' produced no matrix cache; "
        << "'cache path' will be ignored." << std::endl;
  }
  return interpolation;
}

Regridder::Regridder(const eckit::Configuration& conf,
                     const atlas::FunctionSpace& sourceFunctionSpace,
                     const atlas::FunctionSpace& targetFunctionSpace)
    : sourceFunctionSpace_(sourceFunctionSpace),
      targetFunctionSpace_(targetFunctionSpace),
      interpolation_(makeInterpolation(conf, sourceFunctionSpace_,
                                       targetFunctionSpace_)) {}

atlas::FieldSet Regridder::execute(const atlas::FieldSet& source) const {
  atlas::FieldSet target;
  for (atlas::idx_t i = 0; i < source.size(); ++i) {
    target.add(execute(source[i]));
  }
  return target;
}

atlas::Field Regridder::execute(const atlas::Field& source) const {
  atlas::Field target = targetFunctionSpace_.createField(
      atlas::option::name(source.name()) |
      atlas::option::levels(source.levels()) |
      atlas::option::datatype(source.datatype()));

  // Preserve metadata (e.g. missing_value) for downstream masking/IO.
  target.metadata() = source.metadata();

  oops::Log::trace() << "orcamodel::Regridder::execute field '"
      << source.name() << "'" << std::endl;
  interpolation_.execute(source, target);
  return target;
}

}  // namespace orcamodel
