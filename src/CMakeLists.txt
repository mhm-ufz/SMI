list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/cmake-modules)

set(LIB_NAME "${PROJECT_NAME}-lib")

# add sources
file(GLOB sources ./*.f90 ./*.F90)
add_library(${LIB_NAME} ${sources})
target_include_directories(${LIB_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

include(CPM)
# download forces v0.2
set (FORCES_NAME forces)
CPMAddPackage(
  NAME            ${FORCES_NAME}
  GIT_REPOSITORY  https://git.ufz.de/chs/forces.git
  GIT_TAG         v0.2.0
  OPTIONS
    "BUILD_TESTING OFF"
)
target_link_libraries(${LIB_NAME} PUBLIC ${FORCES_NAME})

# control debug messages
option(SMIDEBUG "print debug messages" OFF)
if (SMIDEBUG)
  message(STATUS "SMI: print debug messages")
  target_compile_definitions(${LIB_NAME} PRIVATE SMIDEBUG)
else()
  message(STATUS "SMI: print no debug messages")
endif()

# add all compile options (MPI, OpenMP, Lapack, Coverage)
include(compileoptions)
if (CMAKE_WITH_MPI)
  target_compile_definitions(${LIB_NAME} PRIVATE MPI)
  target_link_libraries(${LIB_NAME} PUBLIC MPI::MPI_Fortran)
endif()
if (CMAKE_WITH_OpenMP)
  target_link_libraries(${LIB_NAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()

# set compiling flags for debug and relese version
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set the preprocessor
  target_compile_definitions(${LIB_NAME} PUBLIC "GFORTRAN")
  target_compile_options(${LIB_NAME} PUBLIC
    -ffree-form -ffixed-line-length-132
    $<$<CONFIG:DEBUG>:-fbacktrace -pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized>
    $<$<CONFIG:RELEASE>:-O3 -fno-peel-loops>
    $<$<BOOL:${CMAKE_WITH_COVERAGE}>:-g -fprofile-arcs -ftest-coverage>
  )
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_definitions(${LIB_NAME} PUBLIC "INTEL")
  target_compile_options(${LIB_NAME} PUBLIC
    -nofixed "SHELL:-assume byterecl" "SHELL:-fp-model source" -m64 "SHELL:-assume realloc-lhs"
    # -fstack-protector-all -fstack-security-check were previously on in debug mode, still needed?
    $<$<CONFIG:DEBUG>:-g "SHELL:-warn all" "SHELL:-check all" -debug -traceback -fp-stack-check -O0>
    $<$<CONFIG:RELEASE>:-O3 -qoverride-limits>
  )
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
  target_compile_definitions(${LIB_NAME} PUBLIC "NAG")
  target_compile_options(${LIB_NAME} PUBLIC
    -fpp -colour -unsharedf95 -ideclient
    # "-C=all" is not set, only "-C -C=alias -C=dangling" and "-ieee=full" instead of "-ieee=stop" because
    # this effectively omits the -C=intovf flag which checks for integer overflow
    # we need to exclude that as the random number generator relies on that technique
    # -ieee=full is needed for mo_utils (is_nan, is_finite etc. fails with -ieee=stop)
    $<$<CONFIG:DEBUG>:-gline -g -nan -O0 -C -C=alias -C=dangling -strict95 -ieee=full>
    $<$<CONFIG:RELEASE>:-O4 -ieee=full>
  )
endif()
