/*
 * (C) British Crown Copyright 2024 Met Office
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
    oops::Log::trace() << "orcamodel::ErrorCovariance created" << std::endl;
}


ErrorCovariance::~ErrorCovariance() {
    oops::Log::trace() << "orcamodel::ErrorCovariance destructed" << std::endl;
}


void ErrorCovariance::linearize(const State &,
                                const Geometry & resol) {
    geom_.reset(new Geometry(resol));
    oops::Log::trace() << "orcamodel::ErrorCovariance linearize" << std::endl;
}


void ErrorCovariance::multiply(const Increment & dxin,
                               Increment & dxout) const {
    oops::Log::trace() << "orcamodel::ErrorCovariance multiply start dxin"
                       <<  dxin << std::endl;

    dxout = dxin;
    std::string err_message =
      "orcamodel::ErrorCovariance::multiply option not implemented";
    throw eckit::NotImplemented(err_message, Here());

    oops::Log::trace() << "orcamodel::ErrorCovariance multiply end dxout"
                       << dxout << std::endl;
}


void ErrorCovariance::inverseMultiply(const Increment & dxin,
                                      Increment & dxout) const {
    dxout = dxin;
    std::string err_message =
            "orcamodel::ErrorCovariance::inverseMultiply notimplemented ";
    throw eckit::NotImplemented(err_message, Here());
    oops::Log::trace() << "orcamodel::ErrorCovariance inverseMultiply"
                       << std::endl;
}


void ErrorCovariance::randomize(Increment & dx) const {
    oops::Log::trace() << "orcamodel::ErrorCovariance randomize" << std::endl;
    std::string err_message =
            "orcamodel::ErrorCovariance::randomise not implemented ";
    throw eckit::NotImplemented(err_message, Here());
}


void ErrorCovariance::print(std::ostream & os) const {
    os << "orcamodel::ErrorCovariance::print not implemented";
}



}  // namespace orcamodel
