/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <cmath>
#include <iostream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "orca-jedi/errorcovariance/ErrorCovariance.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"


// -----------------------------------------------------------------------------
namespace orcamodel {
// -----------------------------------------------------------------------------

ErrorCovariance::ErrorCovariance(const Geometry & resol,
                                 const oops::Variables &,
                                 const eckit::Configuration & conf,
                                 const State &,
                                 const State &) :
    geom_(std::make_shared<const Geometry>(resol)),
    time_(conf.getString("date"))
{
    oops::Log::trace() << "ErrorCovariance(UM) created" << std::endl;
}


ErrorCovariance::~ErrorCovariance() {
    oops::Log::trace() << "ErrorCovariance(UM) destructed" << std::endl;
}


void ErrorCovariance::linearize(const State &,
                                const Geometry & resol) {
    geom_.reset(new Geometry(resol));
    oops::Log::trace() << "ErrorCovariance(UM) linearize" << std::endl;
}


void ErrorCovariance::multiply(const Increment & dxin,
                               Increment & dxout) const {
    oops::Log::trace() << "ErrorCovariance(UM) multiply start dxin"
                       <<  dxin << std::endl;

    dxout = dxin;
    std::string err_message = 
      "umjedi::ErrorCovariance::multiply option not implemented";
    throw eckit::NotImplemented(err_message, Here());

    oops::Log::trace() << "ErrorCovariance(UM) multiply end dxout"
                       << dxout << std::endl;
}


void ErrorCovariance::inverseMultiply(const Increment & dxin,
                                      Increment & dxout) const {
    dxout = dxin;
    std::string err_message =
            "umjedi::ErrorCovariance::inverseMultiply notimplemented ";
    throw eckit::NotImplemented(err_message, Here());
    oops::Log::trace() << "ErrorCovariance(UM) inverseMultiply" << std::endl;
}


void ErrorCovariance::randomize(Increment & dx) const {
    oops::Log::trace() << "ErrorCovariance(UM) randomize" << std::endl;
    std::string err_message =
            "umjedi::ErrorCovariance::randomise not implemented ";
    throw eckit::NotImplemented(err_message, Here());
}


void ErrorCovariance::print(std::ostream & os) const {
    os << "ErrorCovariance::print not implemented";
}



}  // namespace orcamodel
