/* 
 * (C) British Crown Copyright 2017-2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0 
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */ 

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "atlas/interpolation.h"
#include "atlas/functionspace.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "orca-jedi/geometry/Geometry.h"
#include "orca-jedi/state/State.h"

#include "orca-jedi/getvalues/GetValues.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace orcamodel {
  class State;
  class Geometry;

  atlas::functionspace::PointCloud atlasObsFuncSpaceFactory (const ufo::Locations & locs) {

      size_t nlocs = locs.size();

      // Setup observation functionspace
      std::vector<atlas::PointXY> atlasPoints(nlocs);
      oops::Log::trace() << "orcamodel::GetValues:: atlasPoints:" << std::endl;
      for (atlas::idx_t i=0; i < nlocs; ++i ) {
        atlasPoints[i] = atlas::PointXY(locs.lons()[i], locs.lats()[i]);
      }
      return atlas::functionspace::PointCloud(atlasPoints);

  };

  GetValues::GetValues(const Geometry & geom, const ufo::Locations & locs,
            const eckit::Configuration & conf) : 
      atlasObsFuncSpace_(atlasObsFuncSpaceFactory(locs)),
      interpolator_(atlas::option::type("finite-element"),
                    geom.funcSpace(),
                    atlasObsFuncSpace_ ) {

    oops::Log::debug() << "orcamodel::GetValues:: atlasObsFuncSpace_:" << atlasObsFuncSpace_
                       << std::endl;

  };

  void GetValues::fillGeoVaLs(const State& state, const util::DateTime& dt_begin, 
                   const util::DateTime& dt_end, ufo::GeoVaLs& geovals) const {
    std::size_t nlocs = geovals.nlocs();
    std::cout << "statefields print: " << state.stateFields() << std::endl; 

    // dummy "interpolation"
    std::vector<double> vals(nlocs);
    //std::string gv_varname = "sea_ice_category_area_fraction";
    size_t nvars = geovals.getVars().size();
    for ( size_t j=0; j < nvars; ++j) {
      auto gv_varname = geovals.getVars()[j];
      if (!state.variables().has(gv_varname)) {
        std::stringstream err_stream;
        err_stream << "orcajedi::GetValues::fillGeoVals geovals varname \" ";
        err_stream << "\" " << gv_varname << " not found in the model state. " << std::endl;
        err_stream << "    add the variable to the state variables and ";
        err_stream << "add a mapping from the geometry to that variable." << std::endl;
        throw eckit::BadParameter(err_stream.str(), Here());
      }
      
      auto nemo_var_name = state.geometry()->nemo_var_name(gv_varname);
      
      atlas::Field tgt_field = atlasObsFuncSpace_.createField<double>( atlas::option::name( gv_varname ) );
      interpolator_.execute( state.stateFields()[nemo_var_name], tgt_field );
      auto field_view = atlas::array::make_view<double, 1>( tgt_field );
      for (std::size_t i=0; i<nlocs; i++) {
        vals[i] = field_view( i ); 
      }

      geovals.putAtLevel(vals, gv_varname, 0);
      std::vector<double> read_gv(nlocs);
      std::cout << "read_gv size: " << read_gv.size() << std::endl; 
      std::cout << "geovals print: " << geovals << std::endl; 
      geovals.get(read_gv, gv_varname);
      std::cout << "field_view first entry: " << field_view(0) << std::endl; 
      std::cout << "geovals first entry: " << read_gv[0] << std::endl; 
    }
  };

}  // namespace orcamodel




