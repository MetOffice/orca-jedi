
/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <memory>
#include <experimental/filesystem>
#include <netcdf>


#include "eckit/filesystem/PathName.h"
#include "oops/util/ObjectCounter.h"

#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"
#include "atlas/field.h"
#include "atlas/array/ArrayView.h"

namespace orcamodel {

/// \brief Manager of a NEMO field file for reading the pseudo-model state
class NemoFieldReader : private util::ObjectCounter<NemoFieldReader> {
public:
    static const std::string classname() {return "orcamodel::NemoFieldReader";}

    NemoFieldReader( eckit::PathName& filename );

    std::vector<atlas::PointXY> read_locs();
    size_t read_dim_size( const std::string& name );
    std::vector<double> read_surf_var(std::string varname);
    void read_surf_var(std::string varname, 
                       atlas::array::ArrayView<double, 1>& field_view);

private:
    NemoFieldReader() : ncFile() {};
    std::unique_ptr<netCDF::NcFile> ncFile = nullptr; 
};
}  // namespace orcamodel

//------------------------------------------------------------------------------------------------------
