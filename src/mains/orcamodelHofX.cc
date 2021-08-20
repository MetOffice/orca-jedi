/*
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "atlas/library/Library.h"
#include "oops/generic/instantiateModelFactory.h"

#include "oops/runs/HofX4D.h"

#include "ufo/ObsTraits.h"
#include "ufo/instantiateObsFilterFactory.h"

#include "orca-jedi/run/Run.h"
#include "orca-jedi/utilities/OrcaModelTraits.h"

int main(int argc,  char ** argv) {
  orcamodel::Run run(argc, argv);
  oops::instantiateModelFactory<orcamodel::OrcaModelTraits>();
  atlas::Library::instance().initialise();
  ufo::instantiateObsFilterFactory<ufo::ObsTraits>();
  oops::HofX4D<orcamodel::OrcaModelTraits , ufo::ObsTraits> hofx;
  int i = run.execute(hofx);
  atlas::Library::instance().finalise();
  return i;
}
