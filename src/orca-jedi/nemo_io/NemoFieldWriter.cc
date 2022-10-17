/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "NemoFieldWriter.h"

#include <netcdf>
//#include <netcdf> if using Lynton Appel's netcdf-cxx4 from
// https://github.com/Unidata/netcdf-cxx4

#include <algorithm>
#include <sstream>
#include <fstream>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Duration.h"

#include "atlas/runtime/Log.h"
#include "atlas/field.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"

namespace orcamodel {

NemoFieldWriter::NemoFieldWriter(const eckit::PathName& filename, const atlas::Mesh& mesh,
    const std::vector<util::DateTime> datetimes, const std::vector<double> depths)
    : mesh_(mesh)
    , orcaGrid_( mesh.grid() )
    , datetimes_( datetimes )
    , nTimes_( datetimes.size() )
    , depths_( depths )
    , nLevels_( depths.size() )
    {

    bool new_file;
    if (!(filename.exists()) || atlas::mpi::rank() == 0 ) {
      ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(), netCDF::NcFile::replace);
      new_file = true;
    } else {
      ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(), netCDF::NcFile::write);
      new_file = false;
    }
    if (ncFile->isNull()) {
        throw eckit::BadValue("orcamodel::NemoFieldWriter::NemoFieldWriter " + filename + " not found or created");
    }

    if (not orcaGrid_)
       throw eckit::BadValue("orcamodel::NemoFieldWriter:: only writes orca grid data");

    ny_orca_halo_ = orcaGrid_.ny() + orcaGrid_.haloNorth();
    iy_glb_max_ = orcaGrid_.ny() + orcaGrid_.haloNorth() - 1;
    iy_glb_min_ = -orcaGrid_.haloSouth();
    ix_glb_max_ = orcaGrid_.nx() + orcaGrid_.haloEast() - 1;
    ix_glb_min_ = -orcaGrid_.haloWest();
    nx_halo_WE_ = orcaGrid_.nx() + orcaGrid_.haloEast() + orcaGrid_.haloWest();
    ny_halo_NS_ = orcaGrid_.ny() + orcaGrid_.haloNorth() + orcaGrid_.haloSouth();
    glbarray_offset_  = -( nx_halo_WE_ * iy_glb_min_ ) - ix_glb_min_;
    glbarray_jstride_ = nx_halo_WE_;

    if (new_file) {
      setup_dimensions();
    }
    write_dimensions();
}


void NemoFieldWriter::setup_dimensions() {
    size_t nx = nx_halo_WE_;
    size_t ny = ny_halo_NS_;
    auto nxDim = ncFile->addDim("x", nx);
    auto nyDim = ncFile->addDim("y", ny);
    auto nzDim = ncFile->addDim("z", nLevels_);
    auto ntDim = ncFile->addDim("t", nTimes_);
    auto latVar = ncFile->addVar("nav_lat", netCDF::ncDouble, {nyDim, nxDim});
    auto lonVar = ncFile->addVar("nav_lon", netCDF::ncDouble, {nyDim, nxDim});
    auto levVar = ncFile->addVar("z", netCDF::ncDouble, {nzDim});
    auto timeVar = ncFile->addVar("t", netCDF::ncInt, {ntDim});
}

void NemoFieldWriter::write_dimensions() {
    auto lonlat_view = atlas::array::make_view<double,2>(mesh_.nodes().lonlat());
    auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());
    auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));

    auto navLatVar = ncFile->getVar("nav_lat");
    auto navLonVar = ncFile->getVar("nav_lon");
    for (size_t inode = 0; inode < mesh_.nodes().size(); ++inode) {
        const int i = ij(inode, 0) + orcaGrid_.haloWest();
        const int j = ij(inode, 1) + orcaGrid_.haloSouth();
        ATLAS_ASSERT(i >= 0 && i < nx_halo_WE_);
        ATLAS_ASSERT(j >= 0 && j < ny_halo_NS_, "when j is " + std::to_string(j) + " and ny_halo_NS_ " + std::to_string(ny_halo_NS_) );
        std::vector<size_t> f_indices{static_cast<size_t>(j), static_cast<size_t>(i)};
        navLonVar.putVar(f_indices, lonlat_view(inode, 0));
        navLatVar.putVar(f_indices, lonlat_view(inode, 1));
    }

    auto timeVar = ncFile->getVar("t");
    const std::string seconds_since = "seconds since ";
    std::string units_string = "seconds since 1970-01-01 00:00:00";
    timeVar.putAtt("units", units_string);

    units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
    units_string.append("Z");
    util::DateTime epoch = util::DateTime(units_string.substr(seconds_since.size()));
    for (size_t iTime = 0; iTime < nTimes_; ++iTime) {
        int seconds_since_epoch = (datetimes_[iTime] - epoch).toSeconds();
        timeVar.putVar({iTime}, seconds_since_epoch);
    }

    auto levVar = ncFile->getVar("z");
    for (size_t iLevel = 0; iLevel < nLevels_; ++iLevel) {
        levVar.putVar({iLevel}, depths_[iLevel]);
    }
}

void NemoFieldWriter::add_double_variable(std::string& varname) {
    ncFile->addVar(varname, netCDF::ncDouble, (ncFile->getDim("t"), ncFile->getDim("y"), ncFile->getDim("x")));
}

void NemoFieldWriter::write_surf_var(std::string varname, atlas::array::ArrayView<double, 2>& field_view, size_t iTime) {
    try {

        auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());
        auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));

        auto ncVar = ncFile->getVar(varname);
        if (ncVar.isNull()) {
          auto nxDim = ncFile->getDim("x");
          auto nyDim = ncFile->getDim("y");
          auto ntDim = ncFile->getDim("t");
          ncVar = ncFile->addVar(varname, netCDF::ncDouble, {ntDim, nyDim, nxDim});
        }

        auto field_indices = [&](int jLat, int iLon)
        {
          return std::vector<size_t>{static_cast<size_t>(iTime),
                                     static_cast<size_t>(jLat),
                                     static_cast<size_t>(iLon)};
        };

        for (size_t inode = 0; inode<mesh_.nodes().size(); ++inode) {
            //if (ghost(inode)) continue;
            const int i = ij(inode, 0) + orcaGrid_.haloWest();
            const int j = ij(inode, 1) + orcaGrid_.haloSouth();
            ATLAS_ASSERT(i >= 0 && i < nx_halo_WE_);
            ATLAS_ASSERT(j >= 0 && j < ny_halo_NS_);
            ncVar.putVar(field_indices(j, i), field_view(inode, 0));
        }
    }
    catch (netCDF::exceptions::NcException& e) {
        throw eckit::FailedLibraryCall("NetCDF",
            "orcamodel::NemoFieldWriter::write_surf_var", e.what(), Here());
    }
}

void NemoFieldWriter::write_vol_var(std::string varname, atlas::array::ArrayView<double, 2>& field_view, size_t iTime) {
    try {

        auto ghost = atlas::array::make_view<int32_t, 1>(mesh_.nodes().ghost());
        auto ij = atlas::array::make_view<int32_t, 2>(mesh_.nodes().field("ij"));

        auto ncVar = ncFile->getVar(varname);
        if (ncVar.isNull()) {
          auto nxDim = ncFile->getDim("x");
          auto nyDim = ncFile->getDim("y");
          auto nzDim = ncFile->getDim("z");
          auto ntDim = ncFile->getDim("t");
          ncVar = ncFile->addVar(varname, netCDF::ncDouble, {ntDim, nzDim, nyDim, nxDim});
        }

        auto field_indices = [&](int iLev, int jLat, int iLon)
        {
          return std::vector<size_t>{static_cast<size_t>(iTime),
                                     static_cast<size_t>(iLev),
                                     static_cast<size_t>(jLat),
                                     static_cast<size_t>(iLon)};
        };

        for (size_t k = 0; k < nLevels_; ++k) {
            for (size_t inode = 0; inode<mesh_.nodes().size(); ++inode) {
                //if (ghost(inode)) continue;
                const int i = ij(inode, 0) + orcaGrid_.haloWest();
                const int j = ij(inode, 1) + orcaGrid_.haloSouth();
                //std::cout << inode << " ";
                //for (size_t & index : field_indices(j, i)) {std::cout << index << " ";}
                //std::cout << std::endl;
                ncVar.putVar(field_indices(k, j, i), field_view(inode, k));
            }
        }
    }
    catch (netCDF::exceptions::NcException& e) {
        throw eckit::FailedLibraryCall("NetCDF",
            "orcamodel::NemoFieldWriter::write_vol_var", e.what(), Here());
    }
}
}  // namespace orcamodel
