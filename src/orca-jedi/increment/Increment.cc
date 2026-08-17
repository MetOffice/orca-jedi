/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include <string>
#include <vector>
#include <cmath>
#include <ostream>
#include <iostream>
#include <sstream>
#include <limits>
#include <iomanip>

#include "atlas/array/MakeView.h"  // IWYU pragma: keep
#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/field/MissingValue.h"
#include "atlas/functionspace/StructuredColumns.h"  // IWYU pragma: keep

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/log/CodeLocation.h"

#include "oops/base/Variables.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Random.h"

#include "ufo/GeoVaLs.h"

#include "orca-jedi/utilities/IOUtils.h"
#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/increment/IncrementParameters.h"

#include "atlas/mesh.h"  // IWYU pragma: keep
#include "atlas-orca/grid/OrcaGrid.h"

#define INCREMENT_FILL_TOL 1e-6
#define INCREMENT_FILL_VALUE 1e30

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

  std::cout << "Increment(ORCA)::Increment created for " << validTime()
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

/// \brief Copy Increment::Increment constructor.
/// \param other Increment to copy structure from.
/// \param copy Boolean flag copy contents if true.

Increment::Increment(const Increment & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_),
    incrementFields_()
{
  std::cout << "Increment(ORCA)::Increment copy " << copy << std::endl;

  if (copy) {
    incrementFields_ = other.incrementFields_.clone();
  } else {
    incrementFields_ = atlas::FieldSet();
    setupIncrementFields();
  }

  std::cout << "Increment(ORCA)::Increment copied." << std::endl;

  std::cout << "increment copy self print" << std::endl;
  print(std::cout);
  std::cout << "increment copy other print" << std::endl;
  other.print(std::cout);
}

// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------

Increment & Increment::operator=(const Increment & other) {
  std::cout << "Increment(ORCA)::= copy" << std::endl;

  time_ = other.time_;
  vars_ = other.vars_;
  geom_.reset();
  geom_ = other.geom_;
  incrementFields_ = other.incrementFields_.clone();

  std::cout << "Increment(ORCA)::= copy ended" << std::endl;
  return *this;
}

/// += operator: add another increment (dx) to this one
Increment & Increment::operator+=(const Increment & dx) {

  std::cout << "Increment(ORCA)::+= add" << std::endl;

  ASSERT(this->validTime() == dx.validTime());

  /// Marks which nodes are ghost nodes — these must be skipped to avoid double-counting.
  /// ghost is an array of 0 or 1 (for 2 MPI ranks) for every point on this rank's local mesh.
  /// 1 means the point is a ghost node (copy of another rank), 0 means the point j is owned by this rank.

  /// Atlas sets up the ghost flags automatically when the mesh is partitioned across ranks. 
  // The code does not need to know the rank number explicitly — it just trusts ghost(j)

  auto ghost = atlas::array::make_view<int32_t, 1>(geom_->mesh().nodes().ghost());
  
  /// Loop over all fields in the increment and add the corresponding values from dx, skipping ghost nodes and missing values.
  /// Loop over each field pair in the increment and other one (dx)
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field_dx = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName_dx = field_dx.name();
    std::cout << "orcamodel::Increment::add:: field name = " << fieldName
              << " field name dx = " << fieldName_dx
              << std::endl;

    atlas::field::MissingValue mv(field);
    bool has_mv = static_cast<bool>(mv);
    atlas::field::MissingValue mv2(field_dx);
    bool has_mv2 = static_cast<bool>(mv2);

    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view_dx = atlas::array::make_view<double, 2>(field_dx);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          /// Only adds if both values are valid (not missing). If either value is missing, that grid point is left unchanged.
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view_dx(j, k)))) {
              /// the actual addition of the two values at the grid point (j,k)
              field_view(j, k) += field_view_dx(j, k);
              /// this[j,k] = this[j,k] + dx[j,k]
            }
          }
        }
      }
    }
  }

  std::cout << "Increment(ORCA)::+ add ended" << std::endl;
  return *this;
}

/// same for -=
Increment & Increment::operator-=(const Increment & dx) {

  std::cout << "Increment(ORCA)::-= subtract" << std::endl;

  ASSERT(this->validTime() == dx.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field_dx = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName_dx = field_dx.name();
    std::cout << "orcamodel::Increment::subtract:: field name = " << fieldName
              << " field name dx = " << fieldName_dx
              << std::endl;

    atlas::field::MissingValue mv(field);
    bool has_mv = static_cast<bool>(mv);
    atlas::field::MissingValue mv2(field_dx);
    bool has_mv2 = static_cast<bool>(mv2);

    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view_dx = atlas::array::make_view<double, 2>(field_dx);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view_dx(j, k)))) {
              field_view(j, k) -= field_view_dx(j, k);
            }
          }
        }
      }
    }
  }

  std::cout << "Increment(ORCA)::- subtract ended" << std::endl;
  return *this;
}

/// same for *=
Increment & Increment::operator*=(const double & zz) {
  std::cout << "Increment(ORCA)::* multiply" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::cout << "orcamodel::Increment::multiply:: field name = " << fieldName
              << " zz " << zz
              << std::endl;

    atlas::field::MissingValue mv(incrementFields()[fieldName]);
    bool has_mv = static_cast<bool>(mv);

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            field_view(j, k) *= zz;
          }
        }
      }
    }
  }

  std::cout << "Increment(ORCA)::* multiplication ended" << std::endl;
  return *this;
}

/// For owned (valid) points, compute x1 - x2. For ghost points leave both values at 0. 
/// \brief Create increment from the difference of two state objects.
/// \param x1 State object.
/// \param x2 State object subtracted.
void Increment::diff(const State & x1, const State & x2) {

  std::cout << "Increment(ORCA)::diff" << std::endl;

  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field1 = x1.getField(i);
    atlas::Field field2 = x2.getField(i);
    atlas::Field fieldi = incrementFields_[i];

    atlas::field::MissingValue mv1(field1);
    bool has_mv1 = static_cast<bool>(mv1);
    atlas::field::MissingValue mv2(field2);
    bool has_mv2 = static_cast<bool>(mv2);

    std::string fieldName1 = field1.name();
    std::string fieldName2 = field2.name();
    std::string fieldNamei = fieldi.name();
    std::cout << "orcamodel::Increment::diff:: field name 1 = " << fieldName1
              << " field name 2 = " << fieldName2
              << " field name inc = " << fieldNamei
              << std::endl;
    auto field_view1 = atlas::array::make_view<double, 2>(field1);
    auto field_view2 = atlas::array::make_view<double, 2>(field2);
    auto field_viewi = atlas::array::make_view<double, 2>(fieldi);
    for (atlas::idx_t j = 0; j < field_viewi.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_viewi.shape(1); ++k) {
        field_viewi(j, k) = 0;
        if (!ghost(j)) {
          if (!has_mv1 || (has_mv1 && !mv1(field_view1(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view2(j, k)))) {
              field_viewi(j, k) = field_view1(j, k) - field_view2(j, k);
            }
          }
        } else {
          field_viewi(j, k) = 0;
        }
      }
    }
  }
}

/// \brief Set increment fields to a uniform value.
/// \param val Value to use.
void Increment::setval(const double & val) {
  std::cout << "Increment(ORCA)::setval starting" << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::cout << "orcamodel::Increment::setval:: field name = '" << fieldName
              << "' value " << val
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

  std::cout << "Increment(ORCA)::setval done" << std::endl;
}

void Increment::zero() {
  std::cout << "Increment(ORCA)::zero starting" << std::endl;
  this->setval(0);
  std::cout << "Increment(ORCA)::zero done" << std::endl;
}

void Increment::ones() {
  std::cout << "Increment(ORCA)::ones starting" << std::endl;
  this->setval(1);
  std::cout << "Increment(ORCA)::ones done" << std::endl;
}

void Increment::sqrt() {
  std::cout << "Increment(ORCA)::sqrt starting" << std::endl;
  util::sqrtFieldSet(incrementFields_);
  std::cout << "Increment(ORCA)::sqrt done" << std::endl;
}

void Increment::zero(const util::DateTime & vt) {
  time_ = vt;
  std::cout << "orcamodel::Increment::zero at time " << vt << std::endl;
  // NB currently no checking of the time just zeros everything
  this->zero();
}

/// \brief multiply input increment object (x) by a scalar (a) and add onto self (y).
/// \param zz Scalar value (a).
/// \param dx Other increment object (x).
/// \param bool check Check (if true) the validity time of the increments fields matches.
void Increment::axpy(const double & zz, const Increment & dx, const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field_dx = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName_dx = field_dx.name();
    std::cout << "orcamodel::Increment::axpy:: field name = " << fieldName
              << " field name dx = " << fieldName_dx
              << " zz = " << zz
              << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view_dx = atlas::array::make_view<double, 2>(field_dx);

    atlas::field::MissingValue mv(field);
    bool has_mv = static_cast<bool>(mv);
    atlas::field::MissingValue mv2(field_dx);
    bool has_mv2 = static_cast<bool>(mv2);

    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view_dx(j, k)))) {
              field_view(j, k) += zz * field_view_dx(j, k);
            }
          }
        }
      }
    }
  }
}

/// \brief Dot product self increment object with another increment object
/// \param dx Other increment object.
double Increment::dot_product_with(const Increment & dx) const {
  double zz = 0;
  
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  // Deals with multiple fields
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field_dx = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName_dx = field_dx.name();
    std::cout << "orcamodel::Increment::dot_product_with:: field name = " << fieldName
              << " field name dx = " << fieldName_dx
              << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view_dx = atlas::array::make_view<double, 2>(field_dx);

    atlas::field::MissingValue mv(field);
    bool has_mv = static_cast<bool>(mv);
    atlas::field::MissingValue mv2(field_dx);
    bool has_mv2 = static_cast<bool>(mv2);

    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view_dx(j, k)))) {
              zz += field_view(j, k) * field_view_dx(j, k);
            }
          }
        }
      }
    }
  }
  std::cout << "orcamodel::Increment::dot_product_with ended :: zz = " << zz << std::endl;

  return zz;
}

/// \brief Schur product self increment object with another increment object
/// \param dx Other increment object.
void Increment::schur_product_with(const Increment & dx) {
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  for (int i = 0; i< incrementFields_.size(); i++)
  {
    atlas::Field field = incrementFields_[i];
    atlas::Field field_dx = dx.incrementFields_[i];
    std::string fieldName = field.name();
    std::string fieldName_dx = field_dx.name();
    std::cout << "orcamodel::Increment::schur_product_with:: field name = " << fieldName
              << " field name dx = " << fieldName_dx
              << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    auto field_view_dx = atlas::array::make_view<double, 2>(field_dx);

    atlas::field::MissingValue mv(field);
    bool has_mv = static_cast<bool>(mv);
    atlas::field::MissingValue mv2(field_dx);
    bool has_mv2 = static_cast<bool>(mv2);

    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            if (!has_mv2 || (has_mv2 && !mv2(field_view_dx(j, k)))) {
              field_view(j, k) *= field_view_dx(j, k);
            }
          }
        }
      }
    }
  }
}

/// \brief Initialise with a normally distributed random field with a mean of 0 and s.d. of 1.
void Increment::random() {
  std::cout << "orcamodel::Increment::random start" << std::endl;
  std::cout << "orcamodel::Increment::random seed_ " << seed_ << std::endl;

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  
  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::cout << "orcamodel::Increment::random:: field name = " << fieldName
              << std::endl;
    auto field_view = atlas::array::make_view<double, 2>(field);
    // Seed currently hardwired in increment.h
    util::NormalDistribution<double> xx(field_view.shape(0)*field_view.shape(1), 0.0, 1.0, seed_);

    atlas::field::MissingValue mv(incrementFields()[fieldName]);
    bool has_mv = static_cast<bool>(mv);

    int idx = 0;
    for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
      for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
        if (!ghost(j)) {
          if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
            field_view(j, k) = xx[idx];
            idx++;
          }
        }
      }
    }
  }
}

/// \brief Apply Dirac delta functions to configuration specified points.
void Increment::dirac(const eckit::Configuration & config) {
  dirac(oops::validateAndDeserialize<OrcaDiracParameters>(config));
}

/// \brief Apply Dirac delta functions to params specified points.
void Increment::dirac(const OrcaDiracParameters & params) {
  // Adding a delta function at points specified by ixdir, iydir, izdir

  /// Extract x, y, z indices for delta function locations from parameters
  const std::vector<int> & ixdir = params.ixdir;
  const std::vector<int> & iydir = params.iydir;
  const std::vector<int> & izdir = params.izdir;

  ASSERT(ixdir.size() == iydir.size() && ixdir.size() == izdir.size());
  int ndir = ixdir.size();

  /// Get the ORCA grid and compute total width (including halos); prepare storage for flattened node indices.
  atlas::OrcaGrid orcaGrid = geom_->mesh().grid();

  int nx = orcaGrid.nx() + orcaGrid.haloWest() + orcaGrid.haloEast();
  std::vector<int> jpt;

  // Validate each requested point is within field bounds. 
  for (atlas::Field field : incrementFields_) {
    for (int i = 0; i < ndir; i++) {
      if ( (iydir[i]*nx + ixdir[i] >= field.shape(0)) || (izdir[i] >= field.shape(1)) ) {
        std::ostringstream err_stream;
        err_stream << orcamodel::Increment::classname()
                   << " field shape and delta function location configuration mismatch,"
                   << " requested point is out of bounds at: (" << iydir[i]*nx + ixdir[i] << ", "
                   << izdir[i] << ") for field with shape " << field.shape();
        throw eckit::BadValue(err_stream.str(), Here());
      }
    }
  }
  /// Convert 2D (iy, ix) to 1D flattened node index; store and log.
  for (int i = 0; i < ndir; i++) {
    jpt.push_back(iydir[i]*nx + ixdir[i]);
    std::cout << "orcamodel::Increment::dirac:: delta function " << i
              << " at jpt = " << jpt[i]
              << " kpt = " << izdir[i] << std::endl;
  }

  /// Get ghost mask
  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());

  /// zero all fields.
  this->zero();

  /// Loop over fields and requested points; set value to 1 only at owned (non-ghost) nodes, leaving others at 0.

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    std::cout << "orcamodel::Increment::dirac:: field name = " << fieldName
              << std::endl;

    auto field_view = atlas::array::make_view<double, 2>(field);
    for (int i = 0; i < ndir; i++) {
      if (!ghost(jpt[i])) {
        field_view(jpt[i], izdir[i]) = 1;
      }
    }
  }
}


// -----------------------------------------------------------------------------
/// FieldSet operations
// -----------------------------------------------------------------------------

/// \brief Output increment fieldset as an atlas fieldset.
/// \param fset Atlas fieldset to output to.
void Increment::toFieldSet(atlas::FieldSet & fset) const {
  std::cout << "Increment toFieldSet starting" << std::endl;

  fset = atlas::FieldSet();

  for (atlas::idx_t i=0; i < static_cast<atlas::idx_t>(vars_.size()); ++i) {
    // copy variable from increments to new field set
    atlas::Field fieldinc = incrementFields_[i];
    std::string fieldName = fieldinc.name();
    std::cout << "Copy increment toFieldSet " << fieldName << std::endl;

    fset->add(fieldinc);
  }
  std::cout << "Increment toFieldSet done" << std::endl;
}

void Increment::toFieldSetAD(const atlas::FieldSet & fset) {
  std::cout << "Increment toFieldSetAD starting" << std::endl;

  std::string err_message =
      "orcamodel::Increment::toFieldSetAD not implemented";
  throw eckit::NotImplemented(err_message, Here());

  std::cout << "Increment toFieldSetAD done" << std::endl;
}

/// \brief Apply atlas fieldset to an increment fieldset.
/// \param fset Atlas fieldset to apply.
void Increment::fromFieldSet(const atlas::FieldSet & fset) {
  std::cout << "Increment fromFieldSet start" << std::endl;

  for (int i = 0; i< fset.size(); i++) {
    atlas::Field field = fset[i];
    atlas::Field fieldinc = incrementFields_[i];
    std::cout << "Increment fromFieldSet field " << i << " " << field.name() << std::endl;
    std::cout << "Increment fromFieldSet fieldinc " << i
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

  std::cout << "Increment fromFieldSet done" << std::endl;
}

// -----------------------------------------------------------------------------
/// Setup
// -----------------------------------------------------------------------------
/// \brief Setup variables and geometry for increment fields.
void Increment::setupIncrementFields() {
  for (size_t i=0; i < vars_.size(); ++i) {
    // add variable if it isn't already in incrementFields
    std::vector<size_t> varSizes = geom_->variableSizes(vars_);
    if (!incrementFields_.has(vars_[i].name())) {
      atlas::Field field = geom_->functionSpace().createField<double>(
           atlas::option::name(vars_[i].name()) |
           atlas::option::levels(varSizes[i]));
      field.metadata().set("missing_value", INCREMENT_FILL_VALUE);
      field.metadata().set("missing_value_type", "approximately-equals");
      field.metadata().set("missing_value_epsilon", INCREMENT_FILL_TOL);
      field.metadata().set("nearest 3d level", "top");
      incrementFields_.add(field);

      // initialise all data to avoid potential compiler/machine dependent bugs in missingValues
      auto field_view = atlas::array::make_view<double, 2>(field);
      for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
        for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
          field_view(j, k) = 0;
        }
      }

      std::cout << "Increment(ORCA)::setupIncrementFields : "
                << vars_[i].name()
                << " with shape (" << (*(incrementFields_.end()-1)).shape(0)
                << ", " << (*(incrementFields_.end()-1)).shape(1) << ")"
                << std::endl;
      geom_->log_status();
    }
  }
}

// -----------------------------------------------------------------------------
/// I/O
// -----------------------------------------------------------------------------
/// I/O and diagnostics
void Increment::read(const eckit::Configuration & conf) {
  std::string err_message =
      "orcamodel::Increment::read not implemented";
  throw eckit::NotImplemented(err_message, Here());
}

/// \brief Write out increments fields to a file using params specified filename.
void Increment::write(const OrcaIncrementParameters & params) const {
  std::cout << "orcamodel::increment::write" << std::endl;

  std::string output_filename = params.nemoFieldFile.value();
  if (output_filename == "")
    throw eckit::BadValue(std::string("orcamodel::writeIncrementFieldsToFile:: ")
        + "file name not specified", Here());

  auto nemo_field_path = eckit::PathName(output_filename);
  std::cout << "Increment::write to filename "
            << nemo_field_path << std::endl;

  incrementFields_.haloExchange();

  writeFieldsToFile(nemo_field_path, *geom_, time_, incrementFields_);
}

/// \brief Write out increments fields to a file using config specified filename.
void Increment::write(const eckit::Configuration & config) const {
  write(oops::validateAndDeserialize<OrcaIncrementParameters>(config));
}

/// \brief Print some basic information about the self increment object.
void Increment::print(std::ostream & os) const {

  os << "Increment valid at time: " << validTime() << std::endl;
  os << std::string(4, ' ') << vars_ <<  std::endl;
  os << std::string(4, ' ') << "atlas field:" << std::endl;

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    struct Increment::stats s = Increment::stats(fieldName);
    os << std::string(8, ' ') << fieldName <<
          " num: " << s.valid_points <<
          " mean: " << std::setprecision(5) << s.sumx/s.valid_points <<
          " rms: " << std::sqrt(s.sumx2/s.valid_points)  <<
          " min: " << s.min << " max: " << s.max << std::endl;
  }
}

/// \brief Calculate some basic statistics of a field in the increment object.
/// \param fieldName Name of the field to use.
struct Increment::stats Increment::stats(const std::string & fieldName) const {
  struct Increment::stats s;
  s.valid_points = 0;
  s.sumx = 0;
  s.sumx2 = 0;
  s.min = std::numeric_limits<double>::max();
  s.max = std::numeric_limits<double>::lowest();

  auto field_view = atlas::array::make_view<double, 2>(
      incrementFields_[fieldName]);

  auto ghost = atlas::array::make_view<int32_t, 1>(
      geom_->mesh().nodes().ghost());
  atlas::field::MissingValue mv(incrementFields()[fieldName]);
  
  bool has_mv = static_cast<bool>(mv);
  
  for (atlas::idx_t j = 0; j < field_view.shape(0); ++j) {
    for (atlas::idx_t k = 0; k < field_view.shape(1); ++k) {
      if (!ghost(j)) {
        if (!has_mv || (has_mv && !mv(field_view(j, k)))) {
          if (field_view(j, k) > s.max) { s.max=field_view(j, k); }
          if (field_view(j, k) < s.min) { s.min=field_view(j, k); }
          s.sumx += field_view(j, k);
          s.sumx2 += field_view(j, k)*field_view(j, k);
          ++s.valid_points;
        }
      }
    }
  }
  return s;
}

/// \brief Output norm (RMS) of the self increment fields.
double Increment::norm() const {
  int valid_points_all = 0;
  double sumx2all = 0;

  for (atlas::Field field : incrementFields_) {
    std::string fieldName = field.name();
    struct Increment::stats s = Increment::stats(fieldName);
    sumx2all += s.sumx2;
    valid_points_all += s.valid_points;
  }
  // return RMS
  return std::sqrt(sumx2all/valid_points_all);
}

}  // namespace orcamodel
