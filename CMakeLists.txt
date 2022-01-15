cmake_minimum_required(VERSION 3.14)
# additional cmake-modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake-modules)
# get version and date from files (version.txt and version_date.txt)
include(version)
get_version(SMI_VER SMI_VER_DEV SMI_DATE)

# create the project
project(smi
  VERSION ${SMI_VER}
  DESCRIPTION "The Soil Moisture Index - SMI program"
  LANGUAGES Fortran
)
# add full version and date to pre-processor flags (qoutes need in before hand)
add_compile_definitions(PROJVERSION='${SMI_VER_DEV}' PROJDATE='${SMI_DATE}')

# add the executable
add_subdirectory(./app)

# add the library
add_subdirectory(./src)

# setup coverage with GNU
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_WITH_COVERAGE)
  include(CodeCoverage)
  target_link_libraries(smi-lib PUBLIC gcov)
  append_coverage_compiler_flags()
  SETUP_TARGET_FOR_COVERAGE_LCOV(
    NAME smi_coverage_CI
    EXECUTABLE ../CI-scripts/run_cmake_coverage.sh
    DEPENDENCIES smi-exe smi-lib
    GENHTML_ARGS -t "SMI coverage" --html-prolog ../doc/html_files/cov_header.prolog
  )
endif()
