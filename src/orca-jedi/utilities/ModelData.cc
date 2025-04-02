/*
 * (C) Crown Copyright 2024 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <ostream>

#include "orca-jedi/utilities/ModelData.h"


namespace orcamodel {


const eckit::LocalConfiguration ModelData::modelData() const {
  eckit::LocalConfiguration model_data;

  //
  // retrieve data that need to be shared with other system components
  // and stored them into the data structure 'model_data'
  //

  return model_data;
}


void ModelData::print(std::ostream & os) const {
  os << "orcamodel::ModelData::ModelData(): " << modelData();
}


}  // namespace orcamodel

