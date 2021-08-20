/*
 * (C) British Crown Copyright 2017-2018 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/runs/Run.h"

#include "orca-jedi/run/Run.h"

namespace orcamodel {

// -----------------------------------------------------------------------------

Run::Run(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating Run(OrcaModel)" << std::endl;
  const eckit::Configuration * conf = &config();
  //orcamodel_init_f90(&conf);
  oops::Log::trace() << "Run(OrcaModel) created" << std::endl;
}

// -----------------------------------------------------------------------------

Run::~Run() {
  oops::Log::trace() << "Destructing Run(OrcaModel)" << std::endl;
  //orcamodel_finalize_f90();
  //oops::Log::trace() << "Run(OrcaModel): MPI finalized" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace orcamodel
