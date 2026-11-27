
/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <netcdf.h>

#include <string>
#include <vector>

#include "eckit/filesystem/PathName.h"

#include "oops/util/DateTime.h"

namespace orcamodel {

class NemoFieldWriter {
 public:
    static constexpr double fillValue{1.0e+20};
    static const std::string classname() {return "orcamodel::NemoFieldWriter";}

    NemoFieldWriter(const eckit::PathName& filename,
                    const std::vector<util::DateTime>& datetimes,
                    size_t nx, size_t ny,
                    const std::vector<double>& depths);
    ~NemoFieldWriter();
    NemoFieldWriter(const NemoFieldWriter&) = delete;
    NemoFieldWriter& operator=(const NemoFieldWriter&) = delete;
    void write_dimensions(const std::vector<double>& lats,
                          const std::vector<double>& lons);
    template <typename T> void write_surf_var(const std::string& varname,
        const std::vector<T>& var_data, size_t iTime);
    template <typename T> void write_vol_var(const std::string& varname,
        const std::vector<T>& var_data, size_t iTime);

 private:
    void setup_dimensions();
    NemoFieldWriter() = default;
    int ncid_ = -1;
    std::vector<util::DateTime> datetimes_;
    std::vector<double> depths_;
    size_t nLevels_{1};
    size_t nTimes_{1};
    size_t nx_{0};
    size_t ny_{0};
    bool dimension_variables_present_{true};
};
}  // namespace orcamodel

//------------------------------------------------------------------------------------------------------
