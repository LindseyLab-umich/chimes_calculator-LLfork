#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ChimesCalc::ChimesCalc" for configuration "Release"
set_property(TARGET ChimesCalc::ChimesCalc APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ChimesCalc::ChimesCalc PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libchimescalc.so"
  IMPORTED_SONAME_RELEASE "libchimescalc.so"
  )

list(APPEND _cmake_import_check_targets ChimesCalc::ChimesCalc )
list(APPEND _cmake_import_check_files_for_ChimesCalc::ChimesCalc "${_IMPORT_PREFIX}/lib64/libchimescalc.so" )

# Import target "ChimesCalc::ChimesCalc_dynamic" for configuration "Release"
set_property(TARGET ChimesCalc::ChimesCalc_dynamic APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ChimesCalc::ChimesCalc_dynamic PROPERTIES
  IMPORTED_COMMON_LANGUAGE_RUNTIME_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libchimescalc_dl.so"
  IMPORTED_NO_SONAME_RELEASE "TRUE"
  )

list(APPEND _cmake_import_check_targets ChimesCalc::ChimesCalc_dynamic )
list(APPEND _cmake_import_check_files_for_ChimesCalc::ChimesCalc_dynamic "${_IMPORT_PREFIX}/lib64/libchimescalc_dl.so" )

# Import target "ChimesCalc::ChimesCalc_Fortran" for configuration "Release"
set_property(TARGET ChimesCalc::ChimesCalc_Fortran APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ChimesCalc::ChimesCalc_Fortran PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libchimescalc_fortran.so"
  IMPORTED_SONAME_RELEASE "libchimescalc_fortran.so"
  )

list(APPEND _cmake_import_check_targets ChimesCalc::ChimesCalc_Fortran )
list(APPEND _cmake_import_check_files_for_ChimesCalc::ChimesCalc_Fortran "${_IMPORT_PREFIX}/lib64/libchimescalc_fortran.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
