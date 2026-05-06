#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "panvar::panvar" for configuration ""
set_property(TARGET panvar::panvar APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(panvar::panvar PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/bin/panvar"
  )

list(APPEND _cmake_import_check_targets panvar::panvar )
list(APPEND _cmake_import_check_files_for_panvar::panvar "${_IMPORT_PREFIX}/bin/panvar" )

# Import target "panvar::panvarlib" for configuration ""
set_property(TARGET panvar::panvarlib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(panvar::panvarlib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libpanvarlib.a"
  )

list(APPEND _cmake_import_check_targets panvar::panvarlib )
list(APPEND _cmake_import_check_files_for_panvar::panvarlib "${_IMPORT_PREFIX}/lib/libpanvarlib.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
