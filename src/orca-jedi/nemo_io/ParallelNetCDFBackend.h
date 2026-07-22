/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/filesystem/PathName.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"

#include "orca-jedi/nemo_io/IoBackend.h"
#include "orca-jedi/nemo_io/SlabDecomposition.h"

namespace orcamodel {

/// \brief FieldWriteBackend that writes NEMO-layout netCDF files in parallel.
///
/// The file is created collectively on the I/O sub-communicator with
/// nc_create_par (parallel HDF5), and every I/O rank writes its own contiguous
/// y-band hyperslab of each variable using collective netCDF calls. The
/// on-disk layout (dimensions x, y, z, t; coordinate variables nav_lat, nav_lon,
/// t, z; surface variables {t, y, x}; volume variables {t, z, y, x}; fill value
/// NemoFieldWriter::fillValue) is identical to the serial NemoFieldWriter so the
/// two paths are interchangeable.
///
/// All I/O ranks must call every method collectively and in the same order:
/// data variables are defined lazily on first write, which requires a collective
/// define-mode transition across the sub-communicator.
class ParallelNetCDFWriteBackend : public FieldWriteBackend {
 public:
  /// \param io_comm    I/O sub-communicator (from IoPool::io_comm()).
  /// \param path       Output file path.
  /// \param nx, ny     Global grid extents.
  /// \param datetimes  Time coordinate values (defines the t dimension).
  /// \param depths     Depth/level coordinate values (defines the z dimension).
  ParallelNetCDFWriteBackend(const eckit::mpi::Comm& io_comm,
                             const eckit::PathName& path,
                             size_t nx, size_t ny,
                             const std::vector<util::DateTime>& datetimes,
                             const std::vector<double>& depths);

  ~ParallelNetCDFWriteBackend() override;

  ParallelNetCDFWriteBackend(ParallelNetCDFWriteBackend&&) = delete;
  ParallelNetCDFWriteBackend(const ParallelNetCDFWriteBackend&) = delete;
  ParallelNetCDFWriteBackend& operator=(ParallelNetCDFWriteBackend&&) = delete;
  ParallelNetCDFWriteBackend& operator=(const ParallelNetCDFWriteBackend&) = delete;

  void write_dimensions(const Hyperslab& slab,
                        const std::vector<double>& lat_band,
                        const std::vector<double>& lon_band) override;

  void write_surf_slab(const std::string& var_name, size_t t_index,
                       const Hyperslab& slab,
                       const std::vector<double>& data) override;
  void write_surf_slab(const std::string& var_name, size_t t_index,
                       const Hyperslab& slab,
                       const std::vector<float>& data) override;

  void write_vol_slab(const std::string& var_name, size_t t_index,
                      size_t z_begin, size_t z_count, const Hyperslab& slab,
                      const std::vector<double>& data) override;
  void write_vol_slab(const std::string& var_name, size_t t_index,
                      size_t z_begin, size_t z_count, const Hyperslab& slab,
                      const std::vector<float>& data) override;

 private:
  /// \brief Return the varid of \p var_name, defining it collectively (with the
  ///        given netCDF type and dimensionality) on first use.
  int ensure_surf_var(const std::string& var_name, int nc_type);
  int ensure_vol_var(const std::string& var_name, int nc_type);

  int ncid_ = -1;
  int dim_x_ = -1;
  int dim_y_ = -1;
  int dim_z_ = -1;
  int dim_t_ = -1;
  size_t nx_ = 0;
  size_t ny_ = 0;
  size_t n_levels_ = 1;
  size_t n_times_ = 1;
};

/// \brief FieldReadBackend that reads NEMO-layout netCDF files in parallel.
///
/// The mirror image of ParallelNetCDFWriteBackend: the file is opened
/// collectively on the I/O sub-communicator with nc_open_par, and every I/O
/// rank reads its own contiguous y-band hyperslab of a variable using
/// collective netCDF calls. Both 2D surface {t, y, x} and 3D volume
/// {t, z, y, x} variables are supported; the variable's dimensionality is
/// discovered on read.
class ParallelNetCDFReadBackend : public FieldReadBackend {
 public:
  /// \param io_comm I/O sub-communicator (from IoPool::io_comm()).
  /// \param path    Input file path (must be a netCDF-4 / HDF5 file).
  ParallelNetCDFReadBackend(const eckit::mpi::Comm& io_comm,
                            const eckit::PathName& path);

  ~ParallelNetCDFReadBackend() override;

  ParallelNetCDFReadBackend(ParallelNetCDFReadBackend&&) = delete;
  ParallelNetCDFReadBackend(const ParallelNetCDFReadBackend&) = delete;
  ParallelNetCDFReadBackend& operator=(ParallelNetCDFReadBackend&&) = delete;
  ParallelNetCDFReadBackend& operator=(const ParallelNetCDFReadBackend&) = delete;

  void read_slab(const std::string& var_name, size_t t_index, size_t z_index,
                 const Hyperslab& slab, std::vector<double>& data) override;
  void read_slab(const std::string& var_name, size_t t_index, size_t z_index,
                 const Hyperslab& slab, std::vector<float>& data) override;

 private:
  int ncid_ = -1;
};

}  // namespace orcamodel
