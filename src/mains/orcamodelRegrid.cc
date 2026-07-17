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
#include "eckit/filesystem/PathName.h"

#include "oops/runs/Application.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"
#include "orca-jedi/regridder/Regridder.h"
#include "orca-jedi/regridder/SourceExtender.h"
#include "orca-jedi/nemo_io/ReadServer.h"
#include "orca-jedi/utilities/IOUtils.h"
#include "orca-jedi/utilities/MaskUtils.h"

/// \brief Application to regrid ORCA model data to a target grid.
///
/// Reads model state from ORCA grid files, constructs a target grid
/// (e.g. regular lon-lat, rotated pole, or another ORCA grid), and performs
/// interpolation/regridding using the Atlas library.
///
/// An optional source extension step can flood-fill missing-value (land) cells
/// on the source grid before interpolation. This prevents target ocean points
/// from receiving missing values when their entire interpolation stencil falls
/// on source land. The approach is analogous to NEMOVAR's sim_ext module
/// (VAR_SRC/SIM/sim_ext.F90), which iteratively fills coastal land cells with
/// distance-weighted neighbour averages, then smooths the result.
///
/// Configuration is read from a YAML file specifying:
/// - geometry: source ORCA grid configuration
/// - state: source state to regrid
/// - target geometry: target ORCA grid configuration (for ORCA-to-ORCA)
///   OR target grid: atlas grid spec (for ORCA to structured grid)
/// - interpolation method: atlas interpolation configuration
/// - source extension (optional): flood-fill and smoothing parameters
/// - target mask (optional): re-mask regridded fields using a volumetric
///   ancillary field on the target grid (builds a 3D vol_mask)
/// - output: output file path
///
/// Example YAML (ORCA-to-ORCA with source extension and target mask):
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
/// source extension:
///   flood iterations: 2          # cells deep to flood into land (default: 2)
///   smooth iterations: 50        # smoothing passes on flooded cells (default: 50)
///   smooth weight self: 0.35     # self-weight in smoothing (default: 0.35)
///   adjacency: cell-based        # "cell-based" (default) or "edge-based"
/// target mask:                   # optional: re-mask regridded fields on target grid
///   filepath: path/to/target_ancillary.nc
///   variable: votemper           # volumetric field whose missing values define land
/// output:
///   filepath: path/to/output.nc
/// \endcode
///
/// Example YAML (ORCA to structured grid, no source extension):
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
///
/// \note The "source extension" section is optional. Without it, regridding
///   relies solely on Atlas's non_linear missing-value handling, which may
///   leave missing values on target ocean points near coastlines.
///
/// \note The "target mask" section is optional (requires "target geometry").
///   It reads a volumetric field from a NetCDF ancillary on the target grid,
///   derives a 3D land-sea mask (vol_mask) from that field's missing values,
///   and applies it to the regridded result. This ensures that regridded
///   ocean points that should be land on the target grid are correctly set
///   to missing_value at each depth level.
///
/// \note Two adjacency types are supported:
///   - "cell-based" (default): neighbours are all nodes sharing a mesh cell.
///     On a quad mesh this gives 8 neighbours (4 edge + 4 corner), matching
///     the stencil used by NEMOVAR's sim_ext.
///   - "edge-based": neighbours are nodes connected by a mesh edge only.
///     Typically 4–6 neighbours; a sparser but more geometrically principled
///     stencil on unstructured meshes.

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

    // Diagnostics: source fields
    oops::Log::info() << "=== Source field diagnostics ===" << std::endl;
    for (atlas::idx_t i = 0; i < state.stateFields().size(); ++i) {
      const auto& f = state.stateFields()[i];
      oops::Log::info() << "  [" << i << "] name=" << f.name()
                        << " shape=(" << f.shape(0) << "," << f.shape(1) << ")"
                        << " levels=" << f.levels()
                        << " rank=" << f.rank()
                        << " datatype=" << f.datatype().str() << std::endl;
    }

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

    // 3. Optionally extend source fields into land (flood-fill)
    atlas::FieldSet sourceFields;
    for (atlas::idx_t i = 0; i < state.stateFields().size(); ++i) {
      sourceFields.add(state.stateFields()[i].clone());
    }

    if (conf.has("source extension")) {
      const eckit::LocalConfiguration extConf(conf, "source extension");
      int nFlood = extConf.getInt("flood iterations", 2);
      int nSmooth = extConf.getInt("smooth iterations", 50);
      double selfWeight = extConf.getDouble("smooth weight self", 0.35);
      std::string adjStr = extConf.getString("adjacency", "cell-based");
      orcamodel::AdjacencyType adjType =
          (adjStr == "edge-based") ? orcamodel::AdjacencyType::EdgeBased
                                   : orcamodel::AdjacencyType::CellBased;

      orcamodel::SourceExtender extender(
          geom.mesh(), geom.functionSpace(),
          nFlood, nSmooth, selfWeight, adjType);
      geom.log_status();
      extender.extend(sourceFields);
      geom.log_status();

      oops::Log::info() << "Source extension applied: " << nFlood
                        << " flood + " << nSmooth << " smooth iterations ("
                        << adjStr << " adjacency)" << std::endl;
    }

    // 4. Build regridder
    const eckit::LocalConfiguration interpConf(conf, "interpolation method");
    orcamodel::Regridder regridder(interpConf,
                                   geom.functionSpace(),
                                   targetFunctionSpace);
    geom.log_status();

    // 5. Execute regridding
    atlas::FieldSet result = regridder.execute(sourceFields);
    geom.log_status();

    oops::Log::info() << "Regridding complete. Output fields: "
                      << result.size() << std::endl;

    // Diagnostics: regridded result fields
    oops::Log::info() << "=== Regridded result field diagnostics ===" << std::endl;
    for (atlas::idx_t i = 0; i < result.size(); ++i) {
      const auto& f = result[i];
      oops::Log::info() << "  [" << i << "] name=" << f.name()
                        << " shape=(" << f.shape(0) << "," << f.shape(1) << ")"
                        << " levels=" << f.levels()
                        << " rank=" << f.rank()
                        << " datatype=" << f.datatype().str() << std::endl;
    }

    // 5b. Optionally apply target 3D land-sea mask (vol_mask).
    //     Reads a volumetric field from an ancillary on the target grid,
    //     derives a per-level mask from its missing values, and sets
    //     regridded result to missing where the target grid is land.
    if (conf.has("target mask") && targetGeomPtr) {
      const eckit::LocalConfiguration maskConf(conf, "target mask");
      const std::string maskFile = maskConf.getString("filepath");
      const std::string maskVar = maskConf.getString("variable");
      oops::Log::info() << "Applying target mask from '" << maskFile
                        << "' variable '" << maskVar << "'" << std::endl;

      // Read the mask-defining field on the target geometry (all levels)
      const atlas::idx_t nLevels =
          targetGeomPtr->extraFields().field("vol_mask").shape(1);
      atlas::Field maskField = targetGeomPtr->functionSpace().createField<float>(
          atlas::option::name(maskVar) | atlas::option::levels(nLevels));
      {
        orcamodel::ReadServer reader(targetGeomPtr->timer(),
                                     eckit::PathName(maskFile),
                                     targetGeomPtr->mesh());
        auto field_view = atlas::array::make_view<float, 2>(maskField);
        reader.read_var<float>(maskVar, 0, field_view);
        float fill = reader.read_fillvalue<float>(maskVar);
        maskField.metadata().set("missing_value", fill);
        maskField.metadata().set("missing_value_type", "approximately-equals");
        maskField.metadata().set("missing_value_epsilon", 1e-6);
      }

      // Derive vol_mask from the mask field's missing values
      targetGeomPtr->set_vol_mask(maskField);

      // Apply vol_mask to regridded result
      orcamodel::applyMaskToFields(
          targetGeomPtr->extraFields().field("vol_mask"), result);
    }

    // 6. Write output
    if (conf.has("output")) {
      const eckit::LocalConfiguration outputConf(conf, "output");
      const std::string outputPath = outputConf.getString("filepath");

      if (targetGeomPtr) {
        // ORCA target: use existing writeFieldsToFile infrastructure
        const util::DateTime validDate(stateConf.getString("date"));
        eckit::PathName outPath(outputPath);
        if (outPath.exists()) {
          oops::Log::info() << "Removing existing output file: "
                            << outputPath << std::endl;
          outPath.unlink();
        }
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
  // Custom help: intercept --help before oops::Run
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h" || arg == "--help") {
      std::cout <<
R"(OrcaModelRegrid — regrid ORCA model data to a target grid.

Usage:
  orcamodel_regrid.x <config.yaml> [output-log-file]
  orcamodel_regrid.x --help

The YAML configuration file specifies the source geometry, state, target grid,
interpolation method, optional source extension (flood-fill), and output path.

Configuration sections:

  geometry:                        # Source ORCA grid
    grid name: ORCA1_T
    number levels: 3
    nemo variables:
    - {name: sea_water_potential_temperature, nemo field name: votemper,
       model space: volume}

  state:                           # Source state to regrid
    date: 2021-06-30T00:00:00Z
    state variables: [sea_water_potential_temperature]
    nemo field file: /path/to/input.nc

  target geometry:                 # Target ORCA grid (ORCA-to-ORCA)
    grid name: ORCA2_T
    number levels: 3
    nemo variables:
    - {name: sea_water_potential_temperature, nemo field name: votemper,
       model space: volume}

  # OR use 'target grid' for ORCA-to-structured:
  # target grid:
  #   name: L90

  interpolation method:            # Atlas interpolation configuration
    type: unstructured-bilinear-lonlat
    non_linear: missing-if-all-missing

  source extension:                # Optional: flood-fill land cells before regridding
    flood iterations: 2            #   Number of cells to flood into land (default: 2)
    smooth iterations: 50          #   Smoothing passes on flooded cells (default: 50)
    smooth weight self: 0.35       #   Self-weight in smoothing, 0-1 (default: 0.35)
    adjacency: cell-based          #   "cell-based" (default, 8-neighbour, NEMOVAR-like)
                                   #   or "edge-based" (sparser, edge-connected only)

  target mask:                     # Optional: re-mask regridded result on target grid
    filepath: /path/to/target_ancillary.nc  # NetCDF file on the target grid
    variable: votemper             #   Volumetric field whose missing values define land
                                   #   (read at all target levels to build 3D vol_mask)

  output:
    filepath: /path/to/output.nc

Notes:
  - 'target geometry' (ORCA target) or 'target grid' (structured) must be specified.
  - Source extension floods missing-value (land) cells on the source grid with
    neighbour-averaged values before interpolation, preventing missing values on
    target ocean points whose interpolation stencil falls entirely on source land.
    This is analogous to NEMOVAR's sim_ext module.
  - Target mask reads a volumetric field from an ancillary file on the target grid
    and derives a 3D land-sea mask (vol_mask). Where the field has missing values,
    the corresponding (node, level) position is masked. This is applied per-level
    to the regridded result, correctly handling bathymetry differences between grids.
    Requires 'target geometry' (not supported with 'target grid').
  - Output writing is currently only supported for ORCA target grids.
)" << std::endl;
      return 0;
    }
  }

  oops::Run run(argc, argv);
  atlas::Library::instance().initialise();
  OrcaModelRegrid regrid;
  int result = run.execute(regrid);
  atlas::Library::instance().finalise();
  return result;
}
