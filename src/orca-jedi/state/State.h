/*
 * (C) British Crown Copyright 2017-2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>


#include "atlas/field/FieldSet.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/state/StateParameters.h"

namespace ufo {
  class GeoVaLs;
}

namespace oops {
  class Variables;
}

namespace orcamodel {
  class Geometry;
  class GetValuesTraj;
  class Increment;

/// orcaModel model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class State : public util::Printable,
              public util::Serializable,
              private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "orcamodel::State";}

/// Constructor, destructor
  State(const Geometry &,
        const oops::Variables &,
        const util::DateTime &);
  State(const Geometry &,
        const OrcaStateParameters &);
  State(const Geometry &,
        const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const State &);
  virtual ~State();

  State & operator=(const State &);
  void zero();

/// Interactions with Increment
  State & operator+=(const Increment &);

/// I/O and diagnostics
  void read(const OrcaStateParameters &);
  void read(const eckit::Configuration &);
  void analytic_init(const Geometry &);
  void write(const OrcaStateParameters &) const;
  void write(const eckit::Configuration &) const;
  double norm(const std::string & field_name) const;
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

/// Serialize and deserialize
  std::size_t serialSize() const override {return 0;}
  void serialize(std::vector<double> &) const override {}
  void deserialize(const std::vector<double> &, std::size_t &) override {}

/// Utilities
  std::shared_ptr<const Geometry> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

  void updateTime(const util::Duration & dt) {time_ += dt;}


  State(const Geometry & geom, const atlas::FieldSet & fs, util::DateTime & dt)
    : geom_(new Geometry(geom)), vars_(fs.field_names()), time_(dt)
    , stateFields_(fs) {}

  const atlas::FieldSet & stateFields() const {return stateFields_;}
  atlas::FieldSet & stateFields() {return stateFields_;}

  const oops::Variables & variables() const {return vars_;}
  oops::Variables & variables() {return vars_;}

 private:
  void setupStateFields();
  void print(std::ostream &) const override;
  std::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  util::DateTime time_;
  atlas::FieldSet stateFields_;
  OrcaStateParameters params_;
};
// -----------------------------------------------------------------------------

}  // namespace orcamodel
