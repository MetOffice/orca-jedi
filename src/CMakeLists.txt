# (C) British Crown Copyright 2024-2025 Met Office
#

add_subdirectory( orca-jedi )
add_subdirectory( mains )
if( BUILD_TESTING )
  add_subdirectory( tests )
endif()
