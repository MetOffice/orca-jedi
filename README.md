![CI](https://github.com/MetOffice/orca-jedi/actions/workflows/ci.yml/badge.svg)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

&copy; British Crown Copyright 2024 Met Office. All rights reserved.

# orca-jedi

JEDI model interface interface for the NEMO ocean model configurations on ORCA grids from the UK Met Office.

## Description

_orca-jedi_ includes executables to calculate model values at observation locations and to perform quality control. _orca-jedi_ is a JEDI "pseudo-model", meaning that rather than interfacing directly with NEMO or NEMOVAR, the model state is derived from input files. Further applications may be developed based on _orca-jedi_ in the future. (These would be implementations of JEDI OOPS apps, such as for observation generation applications, and various DA applications).

## Getting Started

### Dependencies

  * [cmake](https://cmake.org/)
  * [Unidata/netcdf-cxx4](https://github.com/Unidata/netcdf-cxx4)
  * [ecmwf/ecbuild](https://github.com/ecmwf/ecbuild)
  * [ecmwf/eckit](https://github.com/ecmwf/eckit)
  * [JCSDA/oops](https://github.com/JCSDA/oops)
  * [JCSDA/ufo](https://github.com/JCSDA/ufo)
  * [ecmwf/atlas](https://github.com/ecmwf/atlas)
  * [ecmwf/atlas-orca](https://github.com/ecmwf/atlas-orca)
  * [saber](https://github.com/JCSDA/saber)

### Installing

_orca-jedi_ is a JEDI component package built within [mo-bundle](https://github.com/MetOffice/mo-bundle) at the Met Office - see the README in that project for details on how to build. _orca-jedi_ depends upon the experimental atlas-orca plugin package provided by ECMWF for the atlas library.

Otherwise, it is possible to build _orca-jedi_ within a custom JEDI bundle. For details about JEDI, including installation instructions see the [jedi-docs](http://jedi-docs.jcsda.org/).

These bundles are built, made and installed via cmake, and tested with ctest. All code should be documented at the source level for processing using doxygen. This can be tested by running cmake with the ``-DENABLE_ORCA_JEDI_DOCS=ON`` cmake setting.

### Using the interface

_orca-jedi_ can be run either serially or using MPI. When compiled within the JEDI framework it will produce a series of applications for using with ORCA model data:
```
orcamodelHofX.x <jedi-yaml-config-file>
mpiexec -np 2 orcamodelHofX.x <jedi-yaml-config-file>
```

The jedi configuration is documented in the main jedi documentation. Settings for Met Office operational numerical weather prediction workflows are held internally by the Met Office, however there are some example configurations in the ``examples`` directory as well as the ctest tests inputs (``src/tests/testinputs``).

The parameters for _orca-jedi_ are documented in code in the parameters class. When the program is compiled, these are exported to a yaml schema json file, that can be used in conjunction with your editor to highlight any issues with your configuration. The main sections which likely need to be included regardless of the JEDI application are the ``geometry`` and ``state`` configurations:

```yaml
geometry:
  nemo variables: # list of variables in the geometry
    - name: # internal JEDI name for this variable
      nemo field name: # name in the NEMO field file for the variable
      model space: # "volume", "surface", or "vertical" for 3D, 2D, or 1D field
  grid name: # e.g "eORCA025_T" or "eORCA12_T"
  number levels: # number of vertical levels to include in the grid.
  partitioner: # atlas MPI grid partitioner
state :
  date:  # ISO8601 datetime
  state variables: # list of the model state variables matching the `name` field above
  nemo field file: # input background data
  nemo error field file: # input error data
```

## Help

See the JEDI documentation for help. Additional debugging/trace output is available when:
```
OOPS_DEBUG=true
OOPS_TRACE=true
```

For help with SABER specifically please refer to the [SABER](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/inside/jedi-components/saber/index.html) section in the JEDI documentation. 

## Authors

The current lead maintainer is [@twsearle](https://github.com/twsearle) along with a large amount of help from Met Office contributors (see the "Contributors" page on github).
## Contributing

By contributing you agree to the Contributors License Agreement (CLA) contained in the root directory of the project. Please review this, and if you are able to make a contribution make an issue or pull request for your proposed change. All pull requests should conform to the working practises and be linked to an issue, unless a minor bug fix.

## Working practices

Please see the [JEDI working principles](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/working-practices/index.html) for current working practises.
