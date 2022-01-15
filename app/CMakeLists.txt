# Example application target
add_executable(
  "${PROJECT_NAME}-exe"
  "${CMAKE_CURRENT_SOURCE_DIR}/main.f90"
)
# executable name should be project name
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)
# link against the library in src/
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
# make the app instalable
install(
  TARGETS
  "${PROJECT_NAME}-exe"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
