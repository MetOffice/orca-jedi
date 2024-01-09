
/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <netcdf>

#include <experimental/filesystem>
#include <memory>
#include <string>
#include <vector>

#include "eckit/filesystem/PathName.h"

#include "oops/util/DateTime.h"

#include "atlas/array/ArrayView.h"
#include "atlas/field.h"
#include "atlas/mesh.h"
#include "atlas/runtime/Exception.h"
#include "atlas/util/Point.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

class NemoFieldWriter {
 public:
    NemoFieldWriter(const eckit::PathName& filename, const atlas::Mesh& mesh,
                    const std::vector<util::DateTime> datetimes,
                    const std::vector<double> depths);
    void setup_dimensions();
    void write_dimensions();
    void add_double_variable(std::string& varname);

    void write_surf_var(std::string varname,
        atlas::array::ArrayView<double, 2>& field_view, size_t iTime);
    void write_vol_var(std::string varname,
        atlas::array::ArrayView<double, 2>& field_view, size_t iTime);

 private:
    NemoFieldWriter(): ncFile() {}
    const inline int index_glbarray(int i, int j) {
            ATLAS_ASSERT(i <= ix_glb_max_);
            ATLAS_ASSERT(j <= iy_glb_max_);
            return glbarray_offset_ + j * glbarray_jstride_ + i;
    }

    std::unique_ptr<netCDF::NcFile> ncFile = nullptr;
    atlas::Mesh mesh_;
    atlas::OrcaGrid orcaGrid_;
    std::vector<util::DateTime> datetimes_;
    std::vector<double> depths_;
    size_t nLevels_{1};
    size_t nTimes_{1};
    int ny_orca_halo_{-1};
    int iy_glb_max_{-1};
    int iy_glb_min_{-1};
    int ix_glb_max_{-1};
    int ix_glb_min_{-1};
    int nx_halo_WE_{-1};
    int ny_halo_NS_{-1};
    int glbarray_offset_{-1};
    int glbarray_jstride_{-1};
};
}  // namespace orcamodel

//------------------------------------------------------------------------------------------------------
