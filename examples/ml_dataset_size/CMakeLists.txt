find_package(VTK QUIET)
find_package(Torch QUIET)
if(VTK_FOUND AND Torch_FOUND)
message(STATUS "VTK and Torch were found; enabling support for ml_dataset_size example.")
IRLAddExample_cpp(
  TARGET_NAME
    ml_dataset_size
  HDRS
  SRCS
    "main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/ml_dataset_size"
  COLLECTION
    "examples_ml_classification"
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} ${VTK_CXX_FLAGS}")
target_link_libraries(ml_dataset_size PRIVATE "${TORCH_LIBRARIES}")
target_link_libraries(ml_dataset_size PRIVATE "${VTK_LIBRARIES}")
set_property(TARGET ml_dataset_size PROPERTY CXX_STANDARD 17)
else()
message(STATUS "VTK and/or Torch were not found; disabling support for ml_dataset_size example.")
endif()