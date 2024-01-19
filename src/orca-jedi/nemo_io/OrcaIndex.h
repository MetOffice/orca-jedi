/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <algorithm>
#include <sstream>
#include <limits>
#include <map>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

struct OrcaIndex {
    int32_t ix_glb_max;
    int32_t iy_glb_max;
    int32_t glbarray_offset;
    int32_t glbarray_jstride;
    int32_t nx_halo_WE;
    int32_t ny_halo_NS;

    explicit OrcaIndex(const atlas::OrcaGrid& orcaGrid) {
        iy_glb_max = orcaGrid.ny() + orcaGrid.haloNorth() - 1;
        ix_glb_max = orcaGrid.nx() + orcaGrid.haloEast() - 1;

        nx_halo_WE = orcaGrid.nx() + orcaGrid.haloEast() + orcaGrid.haloWest();
        ny_halo_NS = orcaGrid.ny() + orcaGrid.haloNorth()
          + orcaGrid.haloSouth();

        // vector of local indices: necessary for remote indices of ghost nodes
        int iy_glb_min = -orcaGrid.haloSouth();
        int ix_glb_min = -orcaGrid.haloWest();
        glbarray_offset  = -(nx_halo_WE * iy_glb_min) - ix_glb_min;
        glbarray_jstride = nx_halo_WE;
    }

    int64_t operator()(const int i, const int j) const {
        ATLAS_ASSERT(i <= ix_glb_max,
            std::to_string(i) + " > " + std::to_string(ix_glb_max));
        ATLAS_ASSERT(j <= iy_glb_max,
            std::to_string(j) + " > " + std::to_string(iy_glb_max));
        return glbarray_offset + j * glbarray_jstride + i;
    }
};

}  // namespace orcamodel
