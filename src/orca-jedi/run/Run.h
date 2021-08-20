/*
 * (C) British Crown Copyright 2017-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include "oops/runs/Run.h"

namespace orcamodel {

/*!
 */

// -----------------------------------------------------------------------------

class Run : public oops::Run {
 public:
  Run(int, char **);
  ~Run();
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel
