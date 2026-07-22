/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <netcdf.h>

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "eckit/filesystem/PathName.h"

#include "atlas/parallel/mpi/mpi.h"

#include "oops/util/DateTime.h"

#include "orca-jedi/nemo_io/IoPool.h"
#include "orca-jedi/nemo_io/ParallelNetCDFBackend.h"
#include "orca-jedi/nemo_io/SlabDecomposition.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("ParallelNetCDFWriteBackend writes a NEMO-layout file in parallel") {
  const size_t nx = 10;
  const size_t ny = 7;
  const size_t n_levels = 3;
  const std::vector<double> depths{0.0, 1.0, 2.0};
  const std::vector<util::DateTime> datetimes{util::DateTime("2021-06-30T00:00:00Z")};

  const eckit::PathName path("test_parallel_netcdf.nc");

  // Known-value generators keyed on the flattened global buffer index, so the
  // reader can verify placement independently of the write-side decomposition.
  const auto surf_val = [](size_t g) { return static_cast<double>(g) * 0.5 + 1.0; };
  const auto vol_val = [](size_t k, size_t g) {
    return static_cast<double>(k) * 100.0 + static_cast<double>(g) * 0.25 + 2.0;
  };
  const auto lat_val = [](size_t g) { return static_cast<double>(g) * 0.1; };
  const auto lon_val = [](size_t g) { return static_cast<double>(g) * 0.2; };

  const size_t n_io = std::min<size_t>(atlas::mpi::size(), ny);
  IoPool pool(atlas::mpi::comm(), n_io, nx, ny, "test_pnetcdf_pool");

  {
    if (pool.is_io_rank()) {
      const Hyperslab slab = pool.decomposition().slab(pool.io_rank());
      const size_t g0 = slab.global_index_begin();
      const size_t n = slab.size();

      std::vector<double> lat_band(n), lon_band(n), surf(n);
      for (size_t p = 0; p < n; ++p) {
        lat_band[p] = lat_val(g0 + p);
        lon_band[p] = lon_val(g0 + p);
        surf[p] = surf_val(g0 + p);
      }
      std::vector<double> vol(n_levels * n);
      for (size_t k = 0; k < n_levels; ++k) {
        for (size_t p = 0; p < n; ++p) {
          vol[k * n + p] = vol_val(k, g0 + p);
        }
      }

      ParallelNetCDFWriteBackend backend(pool.io_comm(), path, nx, ny,
                                         datetimes, depths);
      backend.write_dimensions(slab, lat_band, lon_band);
      backend.write_surf_slab("test_surf", 0, slab, surf);
      backend.write_vol_slab("test_vol", 0, 0, n_levels, slab, vol);
    }
  }  // backend destructor closes the file collectively over the I/O comm

  atlas::mpi::comm().barrier();

  // Verify serially on the world root (which is always I/O rank 0).
  if (atlas::mpi::rank() == 0) {
    int ncid = -1;
    EXPECT(nc_open(path.fullName().asString().c_str(), NC_NOWRITE, &ncid) == NC_NOERR);

    int surf_id, vol_id, lat_id, lon_id, z_id;
    EXPECT(nc_inq_varid(ncid, "test_surf", &surf_id) == NC_NOERR);
    EXPECT(nc_inq_varid(ncid, "test_vol", &vol_id) == NC_NOERR);
    EXPECT(nc_inq_varid(ncid, "nav_lat", &lat_id) == NC_NOERR);
    EXPECT(nc_inq_varid(ncid, "nav_lon", &lon_id) == NC_NOERR);
    EXPECT(nc_inq_varid(ncid, "z", &z_id) == NC_NOERR);

    std::vector<double> surf(nx * ny);
    std::vector<double> vol(n_levels * nx * ny);
    std::vector<double> lat(nx * ny), lon(nx * ny), z(n_levels);
    EXPECT(nc_get_var_double(ncid, surf_id, surf.data()) == NC_NOERR);
    EXPECT(nc_get_var_double(ncid, vol_id, vol.data()) == NC_NOERR);
    EXPECT(nc_get_var_double(ncid, lat_id, lat.data()) == NC_NOERR);
    EXPECT(nc_get_var_double(ncid, lon_id, lon.data()) == NC_NOERR);
    EXPECT(nc_get_var_double(ncid, z_id, z.data()) == NC_NOERR);

    for (size_t g = 0; g < nx * ny; ++g) {
      EXPECT(surf[g] == surf_val(g));
      EXPECT(lat[g] == lat_val(g));
      EXPECT(lon[g] == lon_val(g));
    }
    for (size_t k = 0; k < n_levels; ++k) {
      EXPECT(z[k] == depths[k]);
      for (size_t g = 0; g < nx * ny; ++g) {
        EXPECT(vol[k * nx * ny + g] == vol_val(k, g));
      }
    }

    EXPECT(nc_close(ncid) == NC_NOERR);
    std::remove(path.fullName().asString().c_str());
  }
}

//-----------------------------------------------------------------------------

CASE("ParallelNetCDFReadBackend reads NEMO-layout slabs in parallel") {
  const size_t nx = 10;
  const size_t ny = 7;
  const size_t n_levels = 3;
  const std::vector<double> depths{0.0, 1.0, 2.0};
  const std::vector<util::DateTime> datetimes{util::DateTime("2021-06-30T00:00:00Z")};

  const eckit::PathName path("test_parallel_netcdf_read.nc");

  const auto surf_val = [](size_t g) { return static_cast<double>(g) * 0.5 + 1.0; };
  const auto vol_val = [](size_t k, size_t g) {
    return static_cast<double>(k) * 100.0 + static_cast<double>(g) * 0.25 + 2.0;
  };
  const auto lat_val = [](size_t g) { return static_cast<double>(g) * 0.1; };
  const auto lon_val = [](size_t g) { return static_cast<double>(g) * 0.2; };

  const size_t n_io = std::min<size_t>(atlas::mpi::size(), ny);
  IoPool pool(atlas::mpi::comm(), n_io, nx, ny, "test_pnetcdf_read_pool");

  // Write a self-contained netCDF-4 file so the read path never depends on the
  // on-disk format of the pre-existing test data.
  {
    if (pool.is_io_rank()) {
      const Hyperslab slab = pool.decomposition().slab(pool.io_rank());
      const size_t g0 = slab.global_index_begin();
      const size_t n = slab.size();

      std::vector<double> lat_band(n), lon_band(n), surf(n);
      for (size_t p = 0; p < n; ++p) {
        lat_band[p] = lat_val(g0 + p);
        lon_band[p] = lon_val(g0 + p);
        surf[p] = surf_val(g0 + p);
      }
      std::vector<double> vol(n_levels * n);
      for (size_t k = 0; k < n_levels; ++k) {
        for (size_t p = 0; p < n; ++p) {
          vol[k * n + p] = vol_val(k, g0 + p);
        }
      }

      ParallelNetCDFWriteBackend backend(pool.io_comm(), path, nx, ny,
                                         datetimes, depths);
      backend.write_dimensions(slab, lat_band, lon_band);
      backend.write_surf_slab("test_surf", 0, slab, surf);
      backend.write_vol_slab("test_vol", 0, 0, n_levels, slab, vol);
    }
  }

  atlas::mpi::comm().barrier();

  // Read each rank's own slab back through the parallel read backend and check
  // the values against the generators.
  if (pool.is_io_rank()) {
    const Hyperslab slab = pool.decomposition().slab(pool.io_rank());
    const size_t g0 = slab.global_index_begin();
    const size_t n = slab.size();

    ParallelNetCDFReadBackend backend(pool.io_comm(), path);

    std::vector<double> surf;
    backend.read_slab("test_surf", 0, 0, slab, surf);
    EXPECT(surf.size() == n);
    for (size_t p = 0; p < n; ++p) {
      EXPECT(surf[p] == surf_val(g0 + p));
    }

    for (size_t k = 0; k < n_levels; ++k) {
      std::vector<double> lev;
      backend.read_slab("test_vol", 0, k, slab, lev);
      EXPECT(lev.size() == n);
      for (size_t p = 0; p < n; ++p) {
        EXPECT(lev[p] == vol_val(k, g0 + p));
      }
    }
  }

  atlas::mpi::comm().barrier();
  if (atlas::mpi::rank() == 0) {
    std::remove(path.fullName().asString().c_str());
  }
}

//-----------------------------------------------------------------------------

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
  return orcamodel::test::run(argc, argv);
}
