/*
 * (C) British Crown Copyright 2024 MetOffice
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/library/Library.h"
#include "oops/generic/instantiateModelFactory.h"

#include "oops/runs/Run.h"
#include "oops/runs/Variational.h"
#include "saber/oops/instantiateCovarFactory.h"

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
  saber::instantiateCovarFactory<orcamodel::OrcaModelTraits>();
  ufo::instantiateObsFilterFactory();
#if defined(NEMO_FEEDBACK_EXISTS)
  nemo_feedback::instantiateObsFilterFactory();
#endif
  oops::Variational<orcamodel::OrcaModelTraits , ufo::ObsTraits> var;
  int i = run.execute(var);
  atlas::Library::instance().finalise();
  return i;
}
