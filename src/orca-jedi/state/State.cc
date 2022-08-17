/*
 * (C) British Crown Copyright 2017-2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <math.h>

#include <algorithm>
#include <string>
#include <memory>
#include <vector>
#include <functional>
#include <numeric>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/parallel/mpi/mpi.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "eckit/mpi/DataType.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateReadUtils.h"
#include "orca-jedi/model/ModelBias.h"
#include "orca-jedi/model/Model.h"


namespace orcamodel {

// Constructor, destructor
State::State(const Geometry & geom,
             const oops::Variables & vars,
             const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time), params_()
{
  stateFields_ = atlas::FieldSet();

  setupStateFields();

  oops::Log::trace() << "State(ORCA)::State created for "<< validTime()
                     << std::endl;
}

State::State(const Geometry & geom,
             const Parameters_ & params)
  : geom_(new Geometry(geom))
    , vars_(params.stateVariables.value())
    , time_(params.date.value())
    , stateFields_(), params_(params)
{
  std::stringstream params_stream;
  params_stream << "orcamodel::State:: params " << params_;
  oops::Log::debug() << params_stream.str() << std::endl;
  oops::Log::trace() << "State(ORCA)::State:: time: " << validTime()
                     << std::endl;

  setupStateFields();

  if (params_.analyticInit.value().value_or(false)) {
    this->analytic_init(*geom_);
  } else {
    readFieldsFromFile(params_, *geom_, validTime(), "background",
       stateFields_);
    readFieldsFromFile(params_, *geom_, validTime(), "background variance",
       stateFields_);
  }
  oops::Log::trace() << "State(ORCA)::State created." << std::endl;
}

State::State(const Geometry & resol, const State & other)
  : geom_(new Geometry(resol))
    , params_(other.params_)
    , vars_(other.vars_)
    , time_(other.time_)
    , stateFields_(other.stateFields_) {
  ASSERT(other.geom_->grid().uid() == resol.grid().uid());
  oops::Log::trace() << "State(ORCA)::State resolution change: "
                     << " copied as there is no change" << std::endl;
}

State::State(const State & other)
  : geom_(other.geom_)
    , params_(other.params_)
    , vars_(other.vars_)
    , time_(other.time_)
    , stateFields_(other.stateFields_) {
  oops::Log::trace() << "State(ORCA)::State copied." << std::endl;
}

State::~State() {
  oops::Log::trace() << "State(ORCA)::State destructed." << std::endl;
}

// Basic operators

State & State::operator=(const State & rhs) {
  time_ = rhs.time_;
  stateFields_ = rhs.stateFields_;
  vars_ = rhs.vars_;
  geom_.reset();
  geom_ = rhs.geom_;
  return *this;
}

// Interactions with Increments

State & State::operator+=(const Increment & dx) {
  std::string err_message =
      "orcamodel::State::State::operator+= not implemented";
  throw eckit::NotImplemented(err_message, Here());
  oops::Log::trace() << "State(ORCA)::add increment starting" << std::endl;
  oops::Log::trace() << "State(ORCA)::add increment done" << std::endl;
  return *this;
}

// I/O and diagnostics

void State::read(const Parameters_ & params) {
  oops::Log::trace() << "State(ORCA)::read starting for " << params.date.value()
                     << std::endl;

  params_ = params;
  time_ = params.date.value();

  readFieldsFromFile(params, *geom_, validTime(), "background",
      stateFields_);
  oops::Log::trace() << "State(ORCA)::read done" << std::endl;
}

void State::analytic_init(const Geometry & geom) {
  oops::Log::trace() << "State(ORCA)::analytic_init starting" << std::endl;
  this->zero();
  oops::Log::trace() << "State(ORCA)::analytic_init done" << std::endl;
}

void State::setupStateFields() {
  for (size_t i=0; i < vars_.size(); ++i) {
    // add variable if it isn't already in stateFields
    std::vector<size_t> varSizes = geom_->variableSizes(vars_);
    if (!stateFields_.has(vars_[i])) {
      stateFields_.add(geom_->funcSpace().createField<double>(
           atlas::option::name(vars_[i]) |
           atlas::option::levels(varSizes[i])));
    }
  }
}

void State::write(const Parameters_ & params) const {
  oops::Log::trace() << "State(ORCA)::write starting" << std::endl;
  std::string err_message =
      "orcamodel::State::State::write not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

void State::print(std::ostream & os) const {
  oops::Log::trace() << "State(ORCA)::print starting" << std::endl;

  os << std::endl << " Model state valid at time: " << validTime() << std::endl;
  os << std::string(4, ' ') << vars_ <<  std::endl;
  os << std::string(4, ' ') << "atlas field norms:" << std::endl;
  for (atlas::Field field : stateFields_) {
    std::string fieldName = field.name();
    os << std::string(8, ' ') << fieldName << ": " << std::setprecision(5)
       << norm(fieldName) << std::endl;
  }

  oops::Log::trace() << "State(ORCA)::print done" << std::endl;
}

// For accumulator

void State::zero() {
  oops::Log::trace() << "State(ORCA)::zero starting" << std::endl;

  for (atlas::Field field : stateFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::State::zero:: field name = " << fieldName
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        field_view(j, k) = 0;
      }
    }
  }

  oops::Log::trace() << "State(ORCA)::zero done" << std::endl;
}

double State::norm(const std::string & field_name) const {
  double norm = 0;
  int valid_points = 0;

  auto field_view = atlas::array::make_view<double, 2>(
      stateFields_[field_name]);
  auto ghost = atlas::array::make_view<int32_t, 1>(geom_->mesh().nodes().ghost());
  atlas::field::MissingValue mv(stateFields()[field_name]);
  bool has_mv = static_cast<bool>(mv);
  for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
    for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
      if ((!has_mv) || (has_mv && !mv(field_view(j, k)) && !ghost(j))) {
        norm += field_view(j, k)*field_view(j, k);
        ++valid_points;
      }
    }
  }
  return sqrt(norm)/valid_points;
}

}  // namespace orcamodel
