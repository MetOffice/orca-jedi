/*
 * (C) British Crown Copyright 2026 Met Office
 */

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <algorithm>

#include "atlas/field/Field.h"
#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/NodeColumns.h"
#include "atlas/functionspace.h"  // IWYU pragma: keep
#include "atlas/mesh.h"  // IWYU pragma: keep
#include "atlas/grid.h"  // IWYU pragma: keep
#include "atlas/meshgenerator.h"  // IWYU pragma: keep
#include "atlas/parallel/mpi/mpi.h"  // IWYU pragma: keep

#include "eckit/mpi/Comm.h"
#include "eckit/log/Timer.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

#include "orca-jedi/geometry/GeometryParameters.h"
#include "orca-jedi/utilities/Types.h"

namespace atlas {
class Field;
class FieldSet;
class Mesh;
}

namespace orcamodel {

// -----------------------------------------------------------------------------
/// Geometry handles for ORCA model.

  oops::Variables orcaVariableFactory(const eckit::Configuration & config);

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry>{
 public:
  static const std::string classname() {return "orcamodel::Geometry";}

  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  ~Geometry();

  std::vector<size_t> variableSizes(const oops::Variables &) const;
  std::vector<std::string> variableNemoSpaces(const oops::Variables & vars)
      const;
  const eckit::mpi::Comm & getComm() const {return comm_;}
  const oops::Variables & variables() const;
  void create_extrafields();
  void latlon(std::vector<double> & lats, std::vector<double> & lons,
              const bool halo) const;
  const atlas::FunctionSpace & functionSpace() const {return funcSpace_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  const atlas::FieldSet & fields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}

  const atlas::Grid & grid() const {return grid_;}
  const atlas::Mesh & mesh() const {return mesh_;}
  const std::string nemo_var_name(const std::string std_name) const;
  const bool variable_in_variable_type(std::string variable_name,
    std::string variable_type) const;
  bool levelsAreTopDown() const {return true;}
  std::string distributionType() const {
      return params_.partitioner.value();}
  bool parallelOutput() const {return params_.parallelOutput.value();}
  size_t outputIoRanks() const {
      return static_cast<size_t>(std::max(0, params_.outputIoRanks.value()));}
  bool parallelInput() const {return params_.parallelInput.value();}
  size_t inputIoRanks() const {
      return static_cast<size_t>(std::max(0, params_.inputIoRanks.value()));}
  FieldDType fieldPrecision(std::string variable_name) const;
  std::shared_ptr<eckit::Timer> timer() const {return eckit_timer_;}
  void log_status() const;
  /// \brief Record the wall time since the previous phase checkpoint (or since
  ///        the Geometry was constructed) against the named phase and log it to
  ///        oops::Log::info(). Call at the END of a phase; the label names the
  ///        phase that has just completed. Repeated calls with the same label
  ///        accumulate, so read/interp/write totals build up across many field
  ///        reads and writes.
  void log_phase(const std::string & phase) const;
  /// \brief Dump the accumulated per-phase wall times (from log_phase), the
  ///        total elapsed time and the peak resident memory to oops::Log::info().
  void log_phase_summary() const;
  /// \brief As log_phase_summary, but reduce each phase's wall time across the
  ///        supplied communicator with a max (the slowest rank gates collective
  ///        I/O, so the max is the meaningful figure). This call is COLLECTIVE:
  ///        every rank in comm must reach it. Call it at a guaranteed-collective
  ///        point, not from the destructor.
  void log_phase_summary_reduced(const eckit::mpi::Comm & comm) const;
  void set_gmask(atlas::Field &) const;
  void set_vol_mask(atlas::Field &);

 private:
  void print(std::ostream &) const;
  const eckit::mpi::Comm & comm_;
  oops::Variables vars_;
  size_t n_levels_;
  OrcaGeometryParameters params_;
  atlas::Grid grid_;
  atlas::grid::Partitioner partitioner_;
  atlas::Mesh mesh_;
  atlas::functionspace::NodeColumns funcSpace_;
  std::shared_ptr<eckit::Timer> eckit_timer_;
  bool phase_timing_ = false;  ///< Whether per-phase timing/logging is active this run.
  mutable double phase_mark_ = 0.0;  ///< Elapsed time (s) at the last log_phase checkpoint.
  mutable std::map<std::string, double> phase_times_;  ///< Accumulated wall time (s) per phase.
  atlas::FieldSet extraFields_;
};

// -----------------------------------------------------------------------------

}  // namespace orcamodel
