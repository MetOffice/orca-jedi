/*
 * (C) British Crown Copyright 2020 MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/Library.h"
#include "oops/generic/instantiateModelFactory.h"

#include "oops/runs/HofX3D.h"
#include "oops/runs/Run.h"

#include "ufo/ObsTraits.h"
#include "ufo/instantiateObsFilterFactory.h"
#if defined(NEMO_FEEDBACK_EXISTS)
  #include "nemo-feedback/instantiateObsFilterFactory.h"
#endif

#include "orca-jedi/utilities/OrcaModelTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateModelFactory<orcamodel::OrcaModelTraits>();
  atlas::Library::instance().initialise();
  ufo::instantiateObsFilterFactory();
#if defined(NEMO_FEEDBACK_EXISTS)
  nemo_feedback::instantiateObsFilterFactory<ufo::ObsTraits>();
#endif
  oops::HofX3D<orcamodel::OrcaModelTraits , ufo::ObsTraits> hofx;
  int i = run.execute(hofx);
  atlas::Library::instance().finalise();
  return i;
}
