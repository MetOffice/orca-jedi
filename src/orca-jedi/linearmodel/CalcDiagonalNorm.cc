/*
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "orca-jedi/linearmodel/CalcDiagonalNorm.h"

namespace orcamodel {

CalcDiagonalNorm::CalcDiagonalNorm(const orcamodel::State & xx,
                                   const orcamodel::Geometry & tlresol,
                                   const std::string & normMethod,
                                   orcamodel::Increment & invDiagNormMat,
                                   orcamodel::Increment & diagNormMat,
                                   const eckit::LocalConfiguration & conf) {
  oops::Log::error() << "CalcDiagonalNorm::CalcDiagonalNorm not yet implemented" << std::endl;
}

CalcDiagonalNorm::~CalcDiagonalNorm() {}

}  // namespace orcamodel
