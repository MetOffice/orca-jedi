/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "orca-jedi/nemo_io/NemoFieldWriter.h"

// https://github.com/Unidata/netcdf-cxx4
#include <netcdf>

#include <algorithm>
#include <sstream>
#include <fstream>

#include "orca-jedi/utilities/Types.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "atlas/runtime/Log.h"
#include "atlas/field.h"
#include "atlas/array.h"
#include "atlas/parallel/mpi/mpi.h"

namespace orcamodel {

NemoFieldWriter::NemoFieldWriter(const eckit::PathName& filename,
    const std::vector<util::DateTime>& datetimes,
    size_t nx, size_t ny,
    const std::vector<double>& depths) :
      datetimes_(datetimes)
    , nTimes_(datetimes.size())
    , nx_(nx)
    , ny_(ny)
    , depths_(depths)
    , nLevels_(depths.size())
    , dimension_variables_present_(false) {
    bool new_file;
    if (!filename.exists() || atlas::mpi::rank() == 0) {
      ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
          netCDF::NcFile::replace);
      new_file = true;
    } else {
      ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
          netCDF::NcFile::write);
      new_file = false;
    }
    if (ncFile->isNull()) {
        throw eckit::BadValue("orcamodel::NemoFieldWriter::NemoFieldWriter "
            + filename + " not found or created");
    }

    if (new_file) {
      setup_dimensions();
    }

    auto navLatVar = ncFile->getVar("nav_lat");
    auto navLonVar = ncFile->getVar("nav_lon");
    auto tVar = ncFile->getVar("t");
    auto zVar = ncFile->getVar("z");
    if (navLatVar.isNull() || navLonVar.isNull() || tVar.isNull() || zVar.isNull()) {
      dimension_variables_present_ = false;
    } else {
      dimension_variables_present_ = true;
    }
}


void NemoFieldWriter::setup_dimensions() {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::setup_dimensions" << std::endl;
    auto nxDim = ncFile->addDim("x", nx_);
    auto nyDim = ncFile->addDim("y", ny_);
    auto nzDim = ncFile->addDim("z", nLevels_);
    auto ntDim = ncFile->addDim("t", nTimes_);
}

void NemoFieldWriter::write_dimensions(const std::vector<double>& lats,
                                       const std::vector<double>& lons) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_dimensions" << std::endl;

    auto nxDim = ncFile->getDim("x");
    auto nyDim = ncFile->getDim("y");
    auto nzDim = ncFile->getDim("z");
    auto ntDim = ncFile->getDim("t");

    auto navLatVar = ncFile->addVar("nav_lat", netCDF::ncDouble, {nyDim, nxDim});
    auto navLonVar = ncFile->addVar("nav_lon", netCDF::ncDouble, {nyDim, nxDim});
    if ((lons.size() != nx_*ny_) || (lats.size() != nx_*ny_)) {
        throw eckit::BadValue(
            std::string("orcamodel::NemoFieldWriter::write_dimensions")
            + " dimensions sizes do not match lat/lon buffer sizes",
          Here());
    }
    navLonVar.putVar({0, 0}, {ny_, nx_}, lons.data());
    navLatVar.putVar({0, 0}, {ny_, nx_}, lats.data());

    auto timeVar = ncFile->addVar("t", netCDF::ncInt, {ntDim});
    const std::string seconds_since = "seconds since ";
    std::string units_string = "seconds since 1970-01-01 00:00:00";
    timeVar.putAtt("units", units_string);

    units_string.replace(seconds_since.size() + 10, 1, 1, 'T');
    units_string.append("Z");
    util::DateTime epoch = util::DateTime(
        units_string.substr(seconds_since.size()));
    for (size_t iTime = 0; iTime < nTimes_; ++iTime) {
        int seconds_since_epoch = (datetimes_[iTime] - epoch).toSeconds();
        timeVar.putVar({iTime}, seconds_since_epoch);
    }

    auto levVar = ncFile->addVar("z", netCDF::ncDouble, {nzDim});
    for (size_t iLevel = 0; iLevel < nLevels_; ++iLevel) {
        levVar.putVar({iLevel}, depths_[iLevel]);
    }
    dimension_variables_present_ = true;
}
template <typename T> void NemoFieldWriter::write_surf_var(std::string varname,
    const std::vector<T>& var_data, size_t iTime) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_surf_var" << std::endl;
    try {
        if (!dimension_variables_present_) {
          throw eckit::BadValue(
              std::string("orcamodel::NemoFieldWriter::write_surf_var")
              + " can't write '" + varname + "' as the dimensions have not yet been constructed",
            Here());
        }

        auto ncVar = ncFile->getVar(varname);
        if (ncVar.isNull()) {
          auto nxDim = ncFile->getDim("x");
          auto nyDim = ncFile->getDim("y");
          auto ntDim = ncFile->getDim("t");
          ncVar = ncFile->addVar(varname, NetCDFTypeMap<T>::ncType, {ntDim, nyDim, nxDim});
          ncVar.setFill(true, static_cast<T>(NemoFieldWriter::fillValue));
        }

        ncVar.putVar({iTime, 0, 0}, {1, ny_, nx_}, var_data.data());
    }
    catch (netCDF::exceptions::NcException& e) {
        throw eckit::FailedLibraryCall("NetCDF",
            "orcamodel::NemoFieldWriter::write_surf_var", e.what(), Here());
    }
}

template void NemoFieldWriter::write_surf_var<double>(std::string varname,
    const std::vector<double>& var_data, size_t iTime);
template void NemoFieldWriter::write_surf_var<float>(std::string varname,
    const std::vector<float>& var_data, size_t iTime);

template <typename T> void NemoFieldWriter::write_vol_var(std::string varname,
    const std::vector<T>& var_data, size_t iTime) {
    oops::Log::trace() << "orcamodel::NemoFieldWriter::write_vol_var" << std::endl;
    try {
        if (!dimension_variables_present_) {
          throw eckit::BadValue(
              std::string("orcamodel::NemoFieldWriter::write_vol_var")
              + " can't write '" + varname + "' as the dimensions have not yet been constructed",
            Here());
        }

        auto ncVar = ncFile->getVar(varname);
        if (ncVar.isNull()) {
          auto nxDim = ncFile->getDim("x");
          auto nyDim = ncFile->getDim("y");
          auto nzDim = ncFile->getDim("z");
          auto ntDim = ncFile->getDim("t");
          ncVar = ncFile->addVar(varname, NetCDFTypeMap<T>::ncType,
              {ntDim, nzDim, nyDim, nxDim});
          ncVar.setFill(true, static_cast<T>(NemoFieldWriter::fillValue));
        }

        ncVar.putVar({iTime, 0, 0, 0}, {1, nLevels_, ny_, nx_}, var_data.data());
    }
    catch (netCDF::exceptions::NcException& e) {
        throw eckit::FailedLibraryCall("NetCDF",
            "orcamodel::NemoFieldWriter::write_vol_var", e.what(), Here());
    }
}

template void NemoFieldWriter::write_vol_var<double>(std::string varname,
    const std::vector<double>& var_data, size_t iTime);
template void NemoFieldWriter::write_vol_var<float>(std::string varname,
    const std::vector<float>& var_data, size_t iTime);
}  // namespace orcamodel
