/*
 * (C) British Crown Copyright 2026 Met Office
 */

#include "atlas/grid.h"
#include "atlas/functionspace.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/library/Library.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/runs/Application.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/regridder/Regridder.h"
#include "orca-jedi/utilities/IOUtils.h"

/// \brief Application to regrid ORCA model data to a target grid.
///
/// Reads model state from ORCA grid files, constructs a target grid
/// (e.g. regular lon-lat, rotated pole, or another ORCA grid), and performs
/// interpolation/regridding using the Atlas library.
///
/// Configuration is read from a YAML file specifying:
/// - geometry: source ORCA grid configuration
/// - state: source state to regrid
/// - target geometry: target ORCA grid configuration (for ORCA-to-ORCA)
///   OR target grid: atlas grid spec (for ORCA to structured grid)
/// - interpolation method: atlas interpolation configuration
/// - output: output file path
///
/// Example YAML (ORCA-to-ORCA):
/// \code{.yaml}
/// geometry:
///   grid name: ORCA1_T
///   number levels: 3
///   nemo variables:
///   - {name: sea_ice_area_fraction, nemo field name: iiceconc, model space: surface}
///   - {name: sea_water_potential_temperature, nemo field name: votemper, model space: volume}
/// state:
///   date: 2021-06-30T00:00:00Z
///   state variables: [sea_ice_area_fraction, sea_water_potential_temperature]
///   nemo field file: path/to/input.nc
/// target geometry:
///   grid name: ORCA2_T
///   number levels: 3
///   nemo variables:
///   - {name: sea_ice_area_fraction, nemo field name: iiceconc, model space: surface}
///   - {name: sea_water_potential_temperature, nemo field name: votemper, model space: volume}
/// interpolation method:
///   type: unstructured-bilinear-lonlat
///   non_linear: missing-if-all-missing
/// output:
///   filepath: path/to/output.nc
/// \endcode
///
/// Example YAML (ORCA to structured grid):
/// \code{.yaml}
/// geometry:
///   grid name: ORCA2_T
///   number levels: 3
///   nemo variables:
///   - {name: sea_water_potential_temperature, nemo field name: votemper, model space: volume}
/// state:
///   date: 2021-06-30T00:00:00Z
///   state variables: [sea_water_potential_temperature]
///   nemo field file: path/to/input.nc
/// target grid:
///   name: L90
/// interpolation method:
///   type: finite-element
///   non_linear: missing-if-all-missing
/// \endcode

class OrcaModelRegrid : public oops::Application {
 public:
  explicit OrcaModelRegrid(const eckit::mpi::Comm& comm = eckit::mpi::comm())
      : Application(comm) {}

  int execute(const eckit::Configuration& conf) const override {
    oops::Log::info() << "=== OrcaModelRegrid ===" << std::endl;

    // 1. Build source geometry and read state
    const eckit::LocalConfiguration geomConf(conf, "geometry");
    const orcamodel::Geometry geom(geomConf, getComm());

    const eckit::LocalConfiguration stateConf(conf, "state");
    const orcamodel::State state(geom, stateConf);

    oops::Log::info() << "Source state read: " << state << std::endl;

    // 2. Build target function space
    atlas::FunctionSpace targetFunctionSpace;
    std::unique_ptr<orcamodel::Geometry> targetGeomPtr;

    if (conf.has("target geometry")) {
      // ORCA target: build full geometry (enables writing via NemoFieldWriter)
      const eckit::LocalConfiguration targetGeomConf(conf, "target geometry");
      targetGeomPtr = std::make_unique<orcamodel::Geometry>(targetGeomConf, getComm());
      targetFunctionSpace = targetGeomPtr->functionSpace();
      oops::Log::info() << "Target geometry: " << *targetGeomPtr << std::endl;
    } else if (conf.has("target grid")) {
      // Structured grid target
      const eckit::LocalConfiguration targetGridConf(conf, "target grid");
      const std::string targetGridName = targetGridConf.getString("name");
      atlas::Grid targetGrid(targetGridName);
      atlas::MeshGenerator meshgen("structured");
      atlas::Mesh targetMesh = meshgen.generate(targetGrid);
      targetFunctionSpace = atlas::functionspace::NodeColumns(targetMesh);
      oops::Log::info() << "Target grid: " << targetGridName
                        << " (" << targetGrid.size() << " points)" << std::endl;
    } else {
      throw eckit::BadParameter(
          "OrcaModelRegrid: must specify either 'target geometry' or 'target grid'",
          Here());
    }

    // 3. Build regridder
    const eckit::LocalConfiguration interpConf(conf, "interpolation method");
    orcamodel::Regridder regridder(interpConf,
                                   geom.functionSpace(),
                                   targetFunctionSpace);

    // 4. Execute regridding
    atlas::FieldSet result = regridder.execute(state.stateFields());

    oops::Log::info() << "Regridding complete. Output fields: "
                      << result.size() << std::endl;

    // 5. Write output
    if (conf.has("output")) {
      const eckit::LocalConfiguration outputConf(conf, "output");
      const std::string outputPath = outputConf.getString("filepath");

      if (targetGeomPtr) {
        // ORCA target: use existing writeFieldsToFile infrastructure
        const util::DateTime validDate(stateConf.getString("date"));
        orcamodel::writeFieldsToFile(outputPath, *targetGeomPtr, validDate, result);
        oops::Log::info() << "Output written to: " << outputPath << std::endl;
      } else {
        oops::Log::warning() << "Output writing for non-ORCA target grids "
                             << "is not yet implemented. File: " << outputPath
                             << std::endl;
      }
    }

    return 0;
  }

 private:
  std::string appname() const override { return "orcamodel::OrcaModelRegrid"; }
};

int main(int argc, char** argv) {
  oops::Run run(argc, argv);
  atlas::Library::instance().initialise();
  OrcaModelRegrid regrid;
  int result = run.execute(regrid);
  atlas::Library::instance().finalise();
  return result;
}
