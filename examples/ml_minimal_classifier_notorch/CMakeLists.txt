find_package(VTK QUIET)
if(VTK_FOUND)
message(STATUS "VTK was found; enabling support for ml_minimal_classifier_notorch example.")
IRLAddExample_cpp(
  TARGET_NAME
    ml_minimal_classifier_notorch
  HDRS
  SRCS
    "main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/ml_minimal_classifier_notorch"
  COLLECTION
    "examples_ml_minimal_classifier_notorch"
)
target_link_libraries(ml_minimal_classifier_notorch PRIVATE "${VTK_LIBRARIES}")
set_property(TARGET ml_minimal_classifier_notorch PROPERTY CXX_STANDARD 17)
else()
message(STATUS "VTK was not found; disabling support for ml_minimal_classifier_notorch example.")
endif()