/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <ostream>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/log/CodeLocation.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "ufo/GeoVaLs.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/state/StateIOUtils.h"
#include "orca-jedi/increment/Increment.h"

#include "atlas/mesh.h"
#include "atlas-orca/grid/OrcaGrid.h"

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
Increment::Increment(const Geometry & geom,
                     const oops::Variables & vars,
                     const util::DateTime & time)
  : geom_(new Geometry(geom)), vars_(vars), time_(time),
    incrementFields_()
{
  incrementFields_ = atlas::FieldSet();

  setupIncrementFields();

  this->zero();   // may not be needed

  oops::Log::debug() << "Increment(ORCA)::Increment created for "<< validTime()
                     << std::endl;
}

Increment::Increment(const Geometry & geom,
                     const Increment & other)
  : geom_(new Geometry(geom)), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::string err_message =
      "orcamodel::Increment::constructor(geom, other) not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

Increment::Increment(const Increment & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  oops::Log::debug() << "Increment(ORCA)::Increment copy " << copy << std::endl;

  incrementFields_ = atlas::FieldSet();

  setupIncrementFields();

  if (copy) {
    for (size_t i=0; i < vars_.size(); ++i) {
      // copy variable from _Fields to new field set
      atlas::Field field = other.incrementFields_[i];
      oops::Log::debug() << "Copying increment field " << field.name() << std::endl;
      incrementFields_->add(field);
    }
  }

  oops::Log::debug() << "Increment(ORCA)::Increment copied." << std::endl;

  oops::Log::debug() << "increment copy self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment copy other print";
  other.print(oops::Log::debug());
}

// Basic operators
Increment & Increment::operator=(const Increment & rhs) {
  time_ = rhs.time_;
  incrementFields_ = rhs.incrementFields_;
  vars_ = rhs.vars_;
  geom_.reset();
  geom_ = rhs.geom_;

  oops::Log::debug() << "Increment(ORCA)::= copy ended" << std::endl;
  return *this;
}

Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());

  oops::Log::debug() << "increment add self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment add dx print";
  dx.print(oops::Log::debug());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
//  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::add:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) += field_view1(j, k);
      }
    }
  }

  oops::Log::debug() << "increment add self print";
  print(oops::Log::debug());
  oops::Log::debug() << "increment add dx print";
  dx.print(oops::Log::debug());

  oops::Log::debug() << "Increment(ORCA)::+ add ended" << std::endl;
  return *this;
}

Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::subtract:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) -= field_view1(j, k);
      }
    }
  }

  oops::Log::debug() << "Increment(ORCA)::- subtract ended" << std::endl;
  return *this;
}

Increment & Increment::operator*=(const double & zz) {
  oops::Log::debug() << "orcamodel::Increment:multiply start" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::multiply:: field name = " << fieldName
                       << " zz " << zz
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) *= zz;
      }
    }
  }

  oops::Log::debug() << "Increment(ORCA)::* multiplication ended" << std::endl;
  return *this;
}

void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field1 = x1.getField(i);
    atlas::Field field2 = x2.getField(i);
    atlas::Field fieldi = incrementFields_[i];

    std::string fieldName1 = field1.name();
    std::string fieldName2 = field2.name();
    std::string fieldNamei = fieldi.name();
    oops::Log::debug() << "orcamodel::Increment::diff:: field name 1 = " << fieldName1
                       << " field name 2 = " << fieldName2
                       << " field name inc = " << fieldNamei
                       << std::endl;
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    auto field_view2 = atlas::array::make_view<double, 2>(field2);
    auto field_viewi = atlas::array::make_view<double, 2>(fieldi);
    for (atlas::idx_t j = 0; j < field_viewi.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_viewi.shape(1); ++k) {
        if (!ghost(j)) { field_viewi(j, k) = field_view1(j, k) - field_view2(j, k);
        } else { field_viewi(j, k) = 0; }
      }
    }
  }
}

void Increment::setval(const double & val) {
  oops::Log::trace() << "Increment(ORCA)::setval starting" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::setval:: field name = " << fieldName
                       << "value " << val
                       << std::endl;

    atlas::field::MissingValue mv(incrementFields()[fieldName]);
    bool has_mv = static_cast<bool>(mv);

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            field_view(j, k) = val;
          }
        }
      }
    }
  }

  oops::Log::trace() << "Increment(ORCA)::setval done" << std::endl;
}

void Increment::zero() {
  oops::Log::trace() << "Increment(ORCA)::zero starting" << std::endl;
  this->setval(0);
  oops::Log::trace() << "Increment(ORCA)::zero done" << std::endl;
}

void Increment::ones() {
  oops::Log::trace() << "Increment(ORCA)::ones starting" << std::endl;
  this->setval(1);
  oops::Log::trace() << "Increment(ORCA)::ones done" << std::endl;
}

void Increment::zero(const util::DateTime & vt) {
  time_ = vt;
  oops::Log::debug() << "orcamodel::Increment::zero at time " << vt << std::endl;
  // NB currently no checking of the time just zeros everything
  this->zero();
}

void Increment::axpy(const double & zz, const Increment & dx, const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::subtract:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) += zz * field_view1(j, k);
      }
    }
  }
}

double Increment::dot_product_with(const Increment & dx) const {
  double zz = 0;
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  // How should this deal with multiple fields only want to do this with one field??? DJL
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::dot_product_with:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) zz += field_view(j, k) * field_view1(j, k);
      }
    }
  }

  oops::Log::debug() << "orcamodel::Increment::dot_product_with ended :: zz = " << zz << std::endl;

  return zz;
}

void Increment::schur_product_with(const Increment & dx) {
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field1 = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName1 = field1.name();
    oops::Log::debug() << "orcamodel::Increment::schur_product_with:: field name = " << fieldName
                       << " field name 1 = " << fieldName1
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) *= field_view1(j, k);
      }
    }
  }
}

void Increment::random() {
  oops::Log::debug() << "orcamodel::Increment::random start" << std::endl;
  oops::Log::debug() << "orcamodel::Increment::random seed_ " << seed_ << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::random:: field name = " << fieldName
                       << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    // Seed currently hardwired in increment.h DJL - find out how this is dealt with generally
    util::NormalDistribution<double> xx(field_view.shape(0)*field_view.shape(1), 0.0, 1.0, seed_);
    int idx = 0;
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) field_view(j, k) = xx[idx];
        idx++;
      }
    }
  }
}

void Increment::dirac(const eckit::Configuration & conf) {
// Add a delta function at points specified by ixdir, iydir, izdir
  const std::vector<int> & ixdir = conf.getIntVector("ixdir");
  const std::vector<int> & iydir = conf.getIntVector("iydir");
  const std::vector<int> & izdir = conf.getIntVector("izdir");

  ASSERT(ixdir.size() == iydir.size() & ixdir.size() == izdir.size());
  int ndir = ixdir.size();

  atlas::OrcaGrid orcaGrid = geom_->mesh().grid();
  int nx = orcaGrid.nx() + orcaGrid.haloWest() + orcaGrid.haloEast();
  oops::Log::debug() << "orcamodel::Increment::dirac:: nx " << nx << std::endl;

  std::vector<int> jpt;
  for(int i = 0; i < ndir; i++) {
    jpt.push_back(iydir[i]*nx + ixdir[i]);
    oops::Log::debug() << "orcamodel::Increment::dirac:: delta function " << i
                       << " at jpt = " << jpt[i]
                       << " kpt = " << izdir[i] << std::endl;
  }

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  this->zero();

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    oops::Log::debug() << "orcamodel::Increment::dirac:: field name = " << fieldName
                       << std::endl;

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          for (int i = 0; i < ndir; i++) {
            if (j == jpt[i] && k == izdir[i]) {
              field_view(j, k) = 1;
            }
          }
        }
      }
    }
  }
}

void Increment::toFieldSet(atlas::FieldSet & fset) const {
  oops::Log::debug() << "Increment toFieldSet starting" << std::endl;

  fset = atlas::FieldSet();

  for (size_t i=0; i < vars_.size(); ++i) {
    // copy variable from increments to new field set
    atlas::Field fieldinc = incrementFields_[i];
    std::string fieldName = fieldinc.name();
    oops::Log::debug() << "Copy increment toFieldSet " << fieldName << std::endl;

    fset->add(fieldinc);
  }
  oops::Log::debug() << "Increment toFieldSet done" << std::endl;
}

void Increment::toFieldSetAD(const atlas::FieldSet & fset) {
  oops::Log::debug() << "Increment toFieldSetAD starting" << std::endl;

  std::string err_message =
      "orcamodel::Increment::toFieldSetAD not implemented";
  throw eckit::NotImplemented(err_message, Here());

  oops::Log::debug() << "Increment toFieldSetAD done" << std::endl;
}

void Increment::fromFieldSet(const atlas::FieldSet & fset) {
  oops::Log::debug() << "Increment fromFieldSet start" << std::endl;

  for (int i = 0; i< fset.size(); i++) {
    atlas::Field field = fset[i];
    atlas::Field fieldinc = incrementFields_[i];
    oops::Log::debug() << "Increment fromFieldSet field " << i << " " << field.name() << std::endl;
    oops::Log::debug() << "Increment fromFieldSet fieldinc " << i
                       << " " << fieldinc.name() << std::endl;

// copy from field to incrementfields

    auto field_view_to = atlas::array::make_view<double, 2>(fieldinc);
    auto field_view_from = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view_to.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view_to.shape(1); ++k) {
        field_view_to(j, k) = field_view_from(j, k);
      }
    }
  }
  oops::Log::debug() << "Increment fromFieldSet done" << std::endl;
}

void Increment::setupIncrementFields() {
  for (size_t i=0; i < vars_.size(); ++i) {
    // add variable if it isn't already in incrementFields
    std::vector<size_t> varSizes = geom_->variableSizes(vars_);
    if (!incrementFields_.has(vars_[i])) {
      incrementFields_.add(geom_->functionSpace().createField<double>(
           atlas::option::name(vars_[i]) |
           atlas::option::levels(varSizes[i])));
      oops::Log::trace() << "Increment(ORCA)::setupIncrementFields : "
                         << vars_[i] << "has dtype: "
                         << (*(incrementFields_.end()-1)).datatype().str() << std::endl;
      geom_->log_status();
    }
  }
}

/// I/O and diagnostics
void Increment::read(const eckit::Configuration & conf) {
  std::string err_message =
      "orcamodel::Increment::read not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

void Increment::write(const eckit::Configuration & conf) const {
  std::string err_message =
      "orcamodel::Increment::write not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

void Increment::print(std::ostream & os) const {
  double sumx2;
  double sumx;
  double min;
  double max;
  int valid_points;

  oops::Log::trace() << "Increment(ORCA)::print starting" << std::endl;

  os << std::endl << " Increment valid at time: " << validTime() << std::endl;
  os << std::string(4, ' ') << vars_ <<  std::endl;
  os << std::string(4, ' ') << "atlas field:" << std::endl;
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::tie(valid_points, sumx2, sumx, min, max) = stats(fieldName);
    os << std::string(8, ' ') << fieldName <<
          " num: " << valid_points <<
          " mean: " << std::setprecision(5) << sumx/valid_points <<
          " rms: " << sqrt(sumx2/valid_points)  <<
          " min: " << min << " max: " << max << std::endl;
  }
  oops::Log::trace() << "Increment(ORCA)::print done" << std::endl;
}

std::tuple<int, double, double, double, double> Increment::stats(const std::string & fieldName) const {
  int valid_points = 0;
  double sumx = 0;
  double sumx2 = 0;
  double min = 1e30;
  double max = -1e30;

  auto field_view = atlas::array::make_view<double, 2>(
      incrementFields_[fieldName]);
  oops::Log::trace() << "Increment(ORCA):stats" << std::endl;
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  atlas::field::MissingValue mv(incrementFields()[fieldName]);
  bool has_mv = static_cast<bool>(mv);
  for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
    for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
      if (!ghost(j)) {
        if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
          if (field_view(j, k) > max) { max=field_view(j, k); }
          if (field_view(j, k) < min) { min=field_view(j, k); }
          sumx += field_view(j, k);
          sumx2 += field_view(j, k)*field_view(j, k);
          ++valid_points;
        }
      }
    }
  }
  std::cout << "DJL stats " << valid_points << " " << sumx2 << " " << sumx << " " << min << " " << max << std::endl;
  return std::make_tuple(valid_points, sumx2, sumx, min, max);
}

double Increment::norm() const {
  int valid_points;
  int valid_points_all = 0;
  double sumx2all = 0;
  double sumx2;
  double sumx;
  double min;
  double max;

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::tie(valid_points, sumx2, sumx, min, max) = stats(fieldName);
    sumx2all += sumx2;
    valid_points_all += valid_points;
  }
  // return RMS
  return sqrt(sumx2all/valid_points_all);
}

}  // namespace orcamodel
