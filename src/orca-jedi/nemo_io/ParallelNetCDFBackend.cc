/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "orca-jedi/nemo_io/ParallelNetCDFBackend.h"

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>

#include <algorithm>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "orca-jedi/nemo_io/NemoFieldWriter.h"

namespace orcamodel {

namespace {

/// \brief Throw an eckit exception if a netCDF C call did not succeed.
void nc_check(int status, const std::string& where) {
  if (status != NC_NOERR) {
    throw eckit::FailedLibraryCall("NetCDF", where, nc_strerror(status), Here());
  }
}

/// \brief Enable filling with NemoFieldWriter::fillValue for a just-defined
///        variable (must be called in define mode).
void set_fill(int ncid, int varid, int nc_type) {
  if (nc_type == NC_DOUBLE) {
    double fv = static_cast<double>(NemoFieldWriter::fillValue);
    nc_check(nc_def_var_fill(ncid, varid, 0, &fv),
             "ParallelNetCDFWriteBackend::set_fill(double)");
  } else {
    float fv = static_cast<float>(NemoFieldWriter::fillValue);
    nc_check(nc_def_var_fill(ncid, varid, 0, &fv),
             "ParallelNetCDFWriteBackend::set_fill(float)");
  }
}

}  // namespace

ParallelNetCDFWriteBackend::ParallelNetCDFWriteBackend(
    const eckit::mpi::Comm& io_comm, const eckit::PathName& path,
    size_t nx, size_t ny,
    const std::vector<util::DateTime>& datetimes,
    const std::vector<double>& depths)
    : nx_(nx), ny_(ny), n_levels_(depths.size()), n_times_(datetimes.size()),
      n_io_ranks_(io_comm.size()) {
  oops::Log::trace() << "orcamodel::ParallelNetCDFWriteBackend::ctor "
                     << path.fullName().asString() << std::endl;

  // Pick the chunk y-extent to match the y-band decomposition: the number of
  // I/O ranks is clamped to at most ny (a rank cannot own less than one row),
  // and the tallest band is ceil(ny / ranks). Choosing that as the chunk height
  // means each "fat" rank writes exactly one whole chunk row per (t, z) slice,
  // giving contiguous, non-overlapping collective writes that map cleanly onto
  // Lustre stripes. Rows never split, so chunks always span the full x extent.
  const size_t p_eff = std::max<size_t>(1,
      std::min(n_io_ranks_, ny_ == 0 ? size_t{1} : ny_));
  chunk_y_ = std::max<size_t>(1, (ny_ + p_eff - 1) / p_eff);

  // eckit hands out the communicator as a Fortran handle; translate to the C
  // MPI_Comm required by nc_create_par.
  MPI_Comm comm = MPI_Comm_f2c(io_comm.communicator());

  // MPI_Info could carry filesystem striping / collective-buffering hints; the
  // collective HDF5 layer already aggregates writes, so start unhinted.
  nc_check(nc_create_par(path.fullName().asString().c_str(),
                         NC_NETCDF4 | NC_CLOBBER, comm, MPI_INFO_NULL, &ncid_),
           "ParallelNetCDFWriteBackend::nc_create_par");

  // Dimensions, matching NemoFieldWriter::setup_dimensions.
  nc_check(nc_def_dim(ncid_, "x", nx_, &dim_x_), "nc_def_dim x");
  nc_check(nc_def_dim(ncid_, "y", ny_, &dim_y_), "nc_def_dim y");
  nc_check(nc_def_dim(ncid_, "z", n_levels_, &dim_z_), "nc_def_dim z");
  nc_check(nc_def_dim(ncid_, "t", n_times_, &dim_t_), "nc_def_dim t");

  // Coordinate variables (defined now, values written later / below).
  int lat_dims[2] = {dim_y_, dim_x_};
  int lon_dims[2] = {dim_y_, dim_x_};
  int nav_lat_id, nav_lon_id, t_id, z_id;
  nc_check(nc_def_var(ncid_, "nav_lat", NC_DOUBLE, 2, lat_dims, &nav_lat_id),
           "nc_def_var nav_lat");
  nc_check(nc_def_var(ncid_, "nav_lon", NC_DOUBLE, 2, lon_dims, &nav_lon_id),
           "nc_def_var nav_lon");
  nc_check(nc_def_var(ncid_, "t", NC_INT, 1, &dim_t_, &t_id), "nc_def_var t");
  nc_check(nc_def_var(ncid_, "z", NC_DOUBLE, 1, &dim_z_, &z_id), "nc_def_var z");

  // Chunk the 2-D coordinate fields on the same y-band as the data variables.
  size_t coord_chunk[2] = {chunk_y_, nx_};
  nc_check(nc_def_var_chunking(ncid_, nav_lat_id, NC_CHUNKED, coord_chunk),
           "nc_def_var_chunking nav_lat");
  nc_check(nc_def_var_chunking(ncid_, nav_lon_id, NC_CHUNKED, coord_chunk),
           "nc_def_var_chunking nav_lon");

  // Time units attribute (identical string/format to NemoFieldWriter).
  const std::string seconds_since = "seconds since ";
  std::string units_string = "seconds since 1970-01-01 00:00:00";
  nc_check(nc_put_att_text(ncid_, t_id, "units", units_string.size(),
                           units_string.c_str()),
           "nc_put_att_text t units");

  nc_check(nc_enddef(ncid_), "ParallelNetCDFWriteBackend::nc_enddef");

  // Small 1-D coordinate variables carry the same data on every rank, so write
  // them independently from the pool root only.
  nc_check(nc_var_par_access(ncid_, t_id, NC_INDEPENDENT), "par_access t");
  nc_check(nc_var_par_access(ncid_, z_id, NC_INDEPENDENT), "par_access z");
  nc_check(nc_var_par_access(ncid_, nav_lat_id, NC_COLLECTIVE), "par_access nav_lat");
  nc_check(nc_var_par_access(ncid_, nav_lon_id, NC_COLLECTIVE), "par_access nav_lon");

  if (io_comm.rank() == 0) {
    units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
    units_string.append("Z");
    util::DateTime epoch(units_string.substr(seconds_since.size()));
    std::vector<int> t_values(n_times_);
    for (size_t iTime = 0; iTime < n_times_; ++iTime) {
      t_values[iTime] = (datetimes[iTime] - epoch).toSeconds();
    }
    nc_check(nc_put_var_int(ncid_, t_id, t_values.data()), "nc_put_var t");
    nc_check(nc_put_var_double(ncid_, z_id, depths.data()), "nc_put_var z");

    // Report the chunk geometry so it can be matched to Lustre striping. The
    // horizontal footprint (chunk_y * nx) is identical for surface and volume
    // variables because volume levels are written one at a time (z chunk = 1).
    const size_t chunk_elems = chunk_y_ * nx_;
    const size_t f64_bytes = chunk_elems * sizeof(double);
    const size_t f32_bytes = chunk_elems * sizeof(float);
    // Lustre stripe size must be a multiple of 64 KiB; round the f64 chunk up
    // to the next whole MiB for a safe, ready-to-use suggestion.
    const size_t stripe_mib = (f64_bytes + (size_t{1} << 20) - 1) >> 20;
    oops::Log::info()
        << "orcamodel::ParallelNetCDFWriteBackend: parallel output chunking\n"
        << "  output file        : " << path.fullName().asString() << "\n"
        << "  I/O ranks          : " << n_io_ranks_ << "\n"
        << "  global grid (ny,nx): (" << ny_ << ", " << nx_ << ")\n"
        << "  surface chunk shape: {t=1, y=" << chunk_y_ << ", x=" << nx_ << "}\n"
        << "  volume chunk shape : {t=1, z=1, y=" << chunk_y_ << ", x=" << nx_
        << "}\n"
        << "  chunk footprint    : " << chunk_elems << " elements = "
        << f64_bytes << " B (f64) / " << f32_bytes << " B (f32)\n"
        << "  suggested striping : lfs setstripe -c " << n_io_ranks_
        << " -S " << stripe_mib << "M <output-dir>" << std::endl;
  }
}

ParallelNetCDFWriteBackend::~ParallelNetCDFWriteBackend() {
  if (ncid_ >= 0) {
    // Close is collective; ignore the status in the destructor.
    nc_close(ncid_);
  }
}

int ParallelNetCDFWriteBackend::ensure_surf_var(const std::string& var_name,
                                                int nc_type) {
  int varid = -1;
  int status = nc_inq_varid(ncid_, var_name.c_str(), &varid);
  if (status == NC_NOERR) return varid;
  if (status != NC_ENOTVAR) {
    nc_check(status, "ensure_surf_var inq " + var_name);
  }
  int dims[3] = {dim_t_, dim_y_, dim_x_};
  nc_check(nc_redef(ncid_), "ensure_surf_var redef " + var_name);
  nc_check(nc_def_var(ncid_, var_name.c_str(), nc_type, 3, dims, &varid),
           "ensure_surf_var def " + var_name);
  size_t chunk[3] = {1, chunk_y_, nx_};
  nc_check(nc_def_var_chunking(ncid_, varid, NC_CHUNKED, chunk),
           "ensure_surf_var chunking " + var_name);
  set_fill(ncid_, varid, nc_type);
  nc_check(nc_enddef(ncid_), "ensure_surf_var enddef " + var_name);
  nc_check(nc_var_par_access(ncid_, varid, NC_COLLECTIVE),
           "ensure_surf_var par_access " + var_name);
  return varid;
}

int ParallelNetCDFWriteBackend::ensure_vol_var(const std::string& var_name,
                                               int nc_type) {
  int varid = -1;
  int status = nc_inq_varid(ncid_, var_name.c_str(), &varid);
  if (status == NC_NOERR) return varid;
  if (status != NC_ENOTVAR) {
    nc_check(status, "ensure_vol_var inq " + var_name);
  }
  int dims[4] = {dim_t_, dim_z_, dim_y_, dim_x_};
  nc_check(nc_redef(ncid_), "ensure_vol_var redef " + var_name);
  nc_check(nc_def_var(ncid_, var_name.c_str(), nc_type, 4, dims, &varid),
           "ensure_vol_var def " + var_name);
  size_t chunk[4] = {1, 1, chunk_y_, nx_};
  nc_check(nc_def_var_chunking(ncid_, varid, NC_CHUNKED, chunk),
           "ensure_vol_var chunking " + var_name);
  set_fill(ncid_, varid, nc_type);
  nc_check(nc_enddef(ncid_), "ensure_vol_var enddef " + var_name);
  nc_check(nc_var_par_access(ncid_, varid, NC_COLLECTIVE),
           "ensure_vol_var par_access " + var_name);
  return varid;
}

void ParallelNetCDFWriteBackend::write_dimensions(
    const Hyperslab& slab, const std::vector<double>& lat_band,
    const std::vector<double>& lon_band) {
  oops::Log::trace() << "orcamodel::ParallelNetCDFWriteBackend::write_dimensions"
                     << std::endl;
  int nav_lat_id, nav_lon_id;
  nc_check(nc_inq_varid(ncid_, "nav_lat", &nav_lat_id), "inq nav_lat");
  nc_check(nc_inq_varid(ncid_, "nav_lon", &nav_lon_id), "inq nav_lon");

  size_t start[2] = {slab.j_begin, 0};
  size_t count[2] = {slab.j_count, nx_};
  nc_check(nc_put_vara_double(ncid_, nav_lat_id, start, count, lat_band.data()),
           "nc_put_vara nav_lat");
  nc_check(nc_put_vara_double(ncid_, nav_lon_id, start, count, lon_band.data()),
           "nc_put_vara nav_lon");
}

void ParallelNetCDFWriteBackend::write_surf_slab(const std::string& var_name,
    size_t t_index, const Hyperslab& slab, const std::vector<double>& data) {
  const int varid = ensure_surf_var(var_name, NC_DOUBLE);
  size_t start[3] = {t_index, slab.j_begin, 0};
  size_t count[3] = {1, slab.j_count, nx_};
  nc_check(nc_put_vara_double(ncid_, varid, start, count, data.data()),
           "write_surf_slab(double) " + var_name);
}

void ParallelNetCDFWriteBackend::write_surf_slab(const std::string& var_name,
    size_t t_index, const Hyperslab& slab, const std::vector<float>& data) {
  const int varid = ensure_surf_var(var_name, NC_FLOAT);
  size_t start[3] = {t_index, slab.j_begin, 0};
  size_t count[3] = {1, slab.j_count, nx_};
  nc_check(nc_put_vara_float(ncid_, varid, start, count, data.data()),
           "write_surf_slab(float) " + var_name);
}

void ParallelNetCDFWriteBackend::write_vol_slab(const std::string& var_name,
    size_t t_index, size_t z_begin, size_t z_count, const Hyperslab& slab,
    const std::vector<double>& data) {
  const int varid = ensure_vol_var(var_name, NC_DOUBLE);
  size_t start[4] = {t_index, z_begin, slab.j_begin, 0};
  size_t count[4] = {1, z_count, slab.j_count, nx_};
  nc_check(nc_put_vara_double(ncid_, varid, start, count, data.data()),
           "write_vol_slab(double) " + var_name);
}

void ParallelNetCDFWriteBackend::write_vol_slab(const std::string& var_name,
    size_t t_index, size_t z_begin, size_t z_count, const Hyperslab& slab,
    const std::vector<float>& data) {
  const int varid = ensure_vol_var(var_name, NC_FLOAT);
  size_t start[4] = {t_index, z_begin, slab.j_begin, 0};
  size_t count[4] = {1, z_count, slab.j_count, nx_};
  nc_check(nc_put_vara_float(ncid_, varid, start, count, data.data()),
           "write_vol_slab(float) " + var_name);
}

// -----------------------------------------------------------------------------
// Read backend.
// -----------------------------------------------------------------------------

ParallelNetCDFReadBackend::ParallelNetCDFReadBackend(
    const eckit::mpi::Comm& io_comm, const eckit::PathName& path) {
  oops::Log::trace() << "orcamodel::ParallelNetCDFReadBackend::ctor "
                     << path.fullName().asString() << std::endl;
  MPI_Comm comm = MPI_Comm_f2c(io_comm.communicator());
  nc_check(nc_open_par(path.fullName().asString().c_str(), NC_NOWRITE,
                       comm, MPI_INFO_NULL, &ncid_),
           "ParallelNetCDFReadBackend::nc_open_par");
}

ParallelNetCDFReadBackend::~ParallelNetCDFReadBackend() {
  if (ncid_ >= 0) {
    nc_close(ncid_);
  }
}

namespace {

/// \brief Resolve the varid, set collective access and compute the netCDF
///        start/count for a single (t, z) horizontal slice of \p slab. Handles
///        both surface {t, y, x} and volume {t, z, y, x} variables.
int read_slab_setup(int ncid, const std::string& var_name, size_t t_index,
                    size_t z_index, const Hyperslab& slab,
                    size_t start[4], size_t count[4], int& ndims) {
  int varid = -1;
  nc_check(nc_inq_varid(ncid, var_name.c_str(), &varid),
           "ParallelNetCDFReadBackend inq_varid " + var_name);
  nc_check(nc_inq_varndims(ncid, varid, &ndims),
           "ParallelNetCDFReadBackend inq_varndims " + var_name);
  nc_check(nc_var_par_access(ncid, varid, NC_COLLECTIVE),
           "ParallelNetCDFReadBackend par_access " + var_name);
  if (ndims == 3) {
    start[0] = t_index; start[1] = slab.j_begin; start[2] = 0;
    count[0] = 1;       count[1] = slab.j_count; count[2] = slab.nx;
  } else if (ndims == 4) {
    start[0] = t_index; start[1] = z_index; start[2] = slab.j_begin; start[3] = 0;
    count[0] = 1;       count[1] = 1;       count[2] = slab.j_count; count[3] = slab.nx;
  } else {
    throw eckit::FailedLibraryCall("NetCDF",
        "ParallelNetCDFReadBackend::read_slab",
        "variable '" + var_name + "' has unsupported rank", Here());
  }
  return varid;
}

}  // namespace

void ParallelNetCDFReadBackend::read_slab(const std::string& var_name,
    size_t t_index, size_t z_index, const Hyperslab& slab,
    std::vector<double>& data) {
  data.resize(slab.size());
  size_t start[4], count[4];
  int ndims = 0;
  const int varid = read_slab_setup(ncid_, var_name, t_index, z_index, slab,
                                    start, count, ndims);
  nc_check(nc_get_vara_double(ncid_, varid, start, count, data.data()),
           "read_slab(double) " + var_name);
}

void ParallelNetCDFReadBackend::read_slab(const std::string& var_name,
    size_t t_index, size_t z_index, const Hyperslab& slab,
    std::vector<float>& data) {
  data.resize(slab.size());
  size_t start[4], count[4];
  int ndims = 0;
  const int varid = read_slab_setup(ncid_, var_name, t_index, z_index, slab,
                                    start, count, ndims);
  nc_check(nc_get_vara_float(ncid_, varid, start, count, data.data()),
           "read_slab(float) " + var_name);
}

}  // namespace orcamodel
