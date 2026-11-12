/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "orca-jedi/nemo_io/SlabDecomposition.h"

namespace orcamodel {

/// \brief Abstract backend that turns per-I/O-rank slabs of data into bytes in
///        a file.
///
/// This is the seam that decouples orca-jedi's I/O orchestration (communicator
/// pool + redistribution, see IoPool and Redistributor) from the mechanism used
/// to actually write the data. The default implementation is a parallel-netCDF
/// (parallel HDF5) backend in which every I/O rank writes its own contiguous
/// hyperslab collectively. Because the interface only exchanges already
/// redistributed slabs plus their Hyperslab geometry, an alternative backend
/// (for example one delegating to XIOS) can be dropped in without touching the
/// redistribution logic or the callers in WriteServer.
///
/// Only double and float are supported, matching the field element types used
/// throughout orca-jedi; virtual overloads are provided rather than templates
/// so the interface remains polymorphic.
class FieldWriteBackend {
 public:
  virtual ~FieldWriteBackend() = default;

  /// \brief Declare dimensions and write the coordinate (lat/lon) variables.
  /// \param slab      This I/O rank's row-band of the global grid.
  /// \param lat_band  Latitudes for this rank's slab (size == slab.size()).
  /// \param lon_band  Longitudes for this rank's slab (size == slab.size()).
  virtual void write_dimensions(const Hyperslab& slab,
                                const std::vector<double>& lat_band,
                                const std::vector<double>& lon_band) = 0;

  /// \brief Write a 2D surface variable slab at a time index.
  virtual void write_surf_slab(const std::string& var_name, size_t t_index,
                               const Hyperslab& slab,
                               const std::vector<double>& data) = 0;
  virtual void write_surf_slab(const std::string& var_name, size_t t_index,
                               const Hyperslab& slab,
                               const std::vector<float>& data) = 0;

  /// \brief Write a contiguous range of levels of a 3D volume variable slab.
  /// \param z_begin First level in this write.
  /// \param z_count Number of levels in this write (data is laid out level-major:
  ///                data[(k - z_begin) * slab.size() + p]).
  virtual void write_vol_slab(const std::string& var_name, size_t t_index,
                              size_t z_begin, size_t z_count,
                              const Hyperslab& slab,
                              const std::vector<double>& data) = 0;
  virtual void write_vol_slab(const std::string& var_name, size_t t_index,
                              size_t z_begin, size_t z_count,
                              const Hyperslab& slab,
                              const std::vector<float>& data) = 0;
};

/// \brief Abstract backend for reading per-I/O-rank slabs of data from a file.
///
/// The mirror image of FieldWriteBackend: each I/O rank reads its own
/// contiguous hyperslab, which the Redistributor then scatters onto the compute
/// decomposition.
///
/// Reads are expressed as a single horizontal (y, x) slice at a given time and
/// level index, matching the access pattern of ReadServer::read_var (and
/// NemoFieldReader::read_var_slice). The backend inspects the variable's
/// dimensionality: for a 2D surface variable {t, y, x} the level index must be
/// zero and is ignored; for a 3D volume variable {t, z, y, x} it selects the
/// level.
class FieldReadBackend {
 public:
  virtual ~FieldReadBackend() = default;

  /// \brief Read a single horizontal slice of a variable into this rank's slab.
  /// \param var_name Name of the variable in the file.
  /// \param t_index  Time index to read.
  /// \param z_index  Level index (ignored for surface variables).
  /// \param slab     This I/O rank's row-band of the global grid.
  /// \param data     Output slab (resized to slab.size()).
  virtual void read_slab(const std::string& var_name, size_t t_index,
                         size_t z_index, const Hyperslab& slab,
                         std::vector<double>& data) = 0;
  virtual void read_slab(const std::string& var_name, size_t t_index,
                         size_t z_index, const Hyperslab& slab,
                         std::vector<float>& data) = 0;
};

}  // namespace orcamodel
