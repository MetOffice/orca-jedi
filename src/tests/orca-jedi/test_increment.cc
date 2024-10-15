/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <iostream>

#include "eckit/log/Bytes.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

#include "atlas/field/FieldSet.h"
#include "atlas/library/Library.h"

#include "orca-jedi/increment/Increment.h"
#include "orca-jedi/utilities/IOUtils.h"

#include "tests/orca-jedi/OrcaModelTestEnvironment.h"

namespace orcamodel {
namespace test {

//-----------------------------------------------------------------------------

CASE("test increment") {
  EXPECT(eckit::system::Library::exists("atlas-orca"));

  eckit::LocalConfiguration config;
  std::vector<eckit::LocalConfiguration> nemo_var_mappings(4);
  nemo_var_mappings[0].set("name", "sea_ice_area_fraction")
    .set("field precision", "double")
    .set("nemo field name", "iiceconc")
    .set("model space", "surface");
  nemo_var_mappings[1].set("name", "sea_ice_area_fraction_error")
    .set("nemo field name", "sic_tot_var")
    .set("field precision", "double")
    .set("model space", "surface");
  nemo_var_mappings[2].set("name", "sea_surface_foundation_temperature")
    .set("nemo field name", "votemper")
    .set("model space", "surface");
  nemo_var_mappings[3].set("name", "sea_water_potential_temperature")
    .set("nemo field name", "votemper")
    .set("field precision", "double")
    .set("model space", "volume");
  config.set("nemo variables", nemo_var_mappings);
  config.set("grid name", "ORCA2_T");
  config.set("number levels", 10);
  Geometry geometry(config, eckit::mpi::comm());

  eckit::LocalConfiguration inc_config;
  inc_config.set("date", "2021-06-30T00:00:00Z");

  OrcaIncrementParameters params;

  oops::Variables oops_vars2{{oops::Variable{"sea_ice_area_fraction"},
    oops::Variable{"sea_water_potential_temperature"}}};

  oops::Variables oops_vars{{oops::Variable{"sea_ice_area_fraction"}}};

  util::DateTime datetime("2021-06-30T00:00:00Z");

  SECTION("test increment parameters") {
    inc_config.set("output path", "../Data/orca2_t_nemo.nc");
    params.validateAndDeserialize(inc_config);
    EXPECT(params.nemoFieldFile.value() ==
        inc_config.getString("output path"));
    auto datetime = static_cast<util::DateTime>(inc_config.getString("date"));
    EXPECT(params.date.value() == datetime);
  }

  SECTION("test constructor") {
    Increment increment(geometry, oops_vars2, datetime);
    // copy
    Increment increment2 = increment;
  }

  SECTION("test setting increment value") {
    Increment increment(geometry, oops_vars, datetime);
    std::cout << std::endl << "Increment ones: " << std::endl;
    increment.ones();
    increment.print(std::cout);
    EXPECT_EQUAL(increment.norm(), 1);
    std::cout << std::endl << "Increment zero: " << std::endl;
    increment.zero();
    increment.print(std::cout);
    EXPECT_EQUAL(increment.norm(), 0);
    std::cout << std::endl << "Increment random: " << std::endl;
    increment.random();
    increment.print(std::cout);
  }

  SECTION("test dirac") {
    eckit::LocalConfiguration dirac_config;
    std::vector<int> ix = {20, 30};
    std::vector<int> iy = {10, 40};
    std::vector<int> iz = {1, 3};
    dirac_config.set("x indices", ix);
    dirac_config.set("y indices", iy);
    dirac_config.set("z indices", iz);

    Increment increment(geometry, oops_vars, datetime);
    EXPECT_THROWS_AS(increment.dirac(dirac_config), eckit::BadValue);

    dirac_config.set("z indices", std::vector<int>{0, 0});
    increment.dirac(dirac_config);
    increment.print(std::cout);
    EXPECT(std::abs(increment.norm() - 0.0086788) < 1e-6);
  }

  SECTION("test mathematical operators") {
    Increment increment1(geometry, oops_vars, datetime);
    Increment increment2(geometry, oops_vars, datetime);
    increment1.ones();
    std::cout << std::endl << "Increment1:" << std::endl;
    increment1.print(std::cout);
    increment1 *= 3;
    std::cout << std::endl << "Increment1 (*3):" << std::endl;
    increment1.print(std::cout);
    increment2.ones();
    increment2 *= 2;
    std::cout << std::endl << "Increment2 (*2):" << std::endl;
    increment2.print(std::cout);
    increment1 -= increment2;
    std::cout << std::endl << "Increment1 (-increment2):" << std::endl;
    increment1.print(std::cout);
    increment1 += increment2;
    increment1 += increment2;
    std::cout << std::endl << "Increment1 (+increment2*2):" << std::endl;
    increment1.print(std::cout);
    double zz = increment1.dot_product_with(increment2);
    std::cout << std::endl << "Dot product increment1.increment2 = " << zz << std::endl;
    increment1.schur_product_with(increment2);
    std::cout << std::endl << "Increment1 (Schur product with increment 2):" << std::endl;
    increment1.print(std::cout);
    increment1.axpy(100, increment2, true);
    std::cout << std::endl << "Increment1 axpy (increment2*100 + increment1):" << std::endl;
    increment1.print(std::cout);
    EXPECT_EQUAL(increment1.norm(), 210);
  }

  SECTION("test increments to fieldset and back to increments") {
    Increment increment1(geometry, oops_vars, datetime);
    increment1.ones();
    Increment increment2(geometry, oops_vars, datetime);
    increment2.zero();
    atlas::FieldSet incfset = atlas::FieldSet();
    increment1.Increment::toFieldSet(incfset);
    increment2.Increment::fromFieldSet(incfset);
    increment1.print(std::cout);
    increment2.print(std::cout);
    EXPECT_EQUAL(increment1.norm(), increment2.norm());
  }

  SECTION("test increment diff with state inputs") {
    // Using the same variables and double type as the increments
    // Code to deal with differing variables in state and increment not currently implemented
    orcamodel::State state1(geometry, oops_vars, datetime);
    orcamodel::State state2(geometry, oops_vars, datetime);
    state1.zero();
    state2.zero();
    std::cout << "state1 norm:" << oops_vars[0].name();
    std::cout << state1.norm<double>(oops_vars[0].name()) << std::endl;
    std::cout << "state2 norm:" << oops_vars[0].name();
    std::cout << state2.norm<double>(oops_vars[0].name()) << std::endl;
    Increment increment(geometry, oops_vars, datetime);
    increment.diff(state1, state2);
    std::cout << "increment (diff state1 state2):" << std::endl;
    increment.print(std::cout);
  }

  SECTION("test increment write") {
    Increment increment1(geometry, oops_vars, datetime);
    increment1.ones();
    inc_config.set("output path", "../testoutput/orca2_t_increment_output.nc");
    params.validateAndDeserialize(inc_config);
    increment1.write(params);
  }

  SECTION("test orca atlas fieldset write") {
    Increment increment1(geometry, oops_vars, datetime);
    increment1.ones();
    increment1 *=2;
    atlas::FieldSet incfset = atlas::FieldSet();
    increment1.Increment::toFieldSet(incfset);
    writeFieldsToFile("../testoutput/orca2_t_increment_fieldset_output.nc",
        geometry, datetime, incfset);
  }
}

}  // namespace test
}  // namespace orcamodel

int main(int argc, char** argv) {
    return orcamodel::test::run(argc, argv);
}
