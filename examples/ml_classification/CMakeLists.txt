find_package(VTK QUIET)
find_package(Torch QUIET)
if(VTK_FOUND AND Torch_FOUND)
message(STATUS "VTK and Torch were found; enabling support for ml_classification example.")
IRLAddExample_cpp(
  TARGET_NAME
    ml_classification
  HDRS
  SRCS
    "main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/ml_classification"
  COLLECTION
    "examples_ml_classification"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} ${VTK_CXX_FLAGS}")
target_link_libraries(ml_classification PRIVATE "${TORCH_LIBRARIES}")
target_link_libraries(ml_classification PRIVATE "${VTK_LIBRARIES}")
set_property(TARGET ml_classification PROPERTY CXX_STANDARD 17)
else()
message(STATUS "VTK and/or Torch were not found; disabling support for ml_classification example.")
endif()