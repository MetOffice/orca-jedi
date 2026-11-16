/*
 * (C) Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>

#include "oops/base/ModelBase.h"
#include "oops/interface/Model.h"

#include "orca-jedi/utilities/OrcaModelTraits.h"
#include "orca-jedi/model/Model.h"

namespace orcamodel {

void OrcaModelTraits::registerModels() {
  static oops::ModelMaker<OrcaModelTraits,
                          oops::interface::Model<OrcaModelTraits, orcamodel::Model> >
    makerLfricjedi_("OrcaModel");
}

}  // namespace orcamodel
