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
set (Torch_DIR "/home/mlcfd/Repositories/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
target_link_libraries(ml_classification PRIVATE "${TORCH_LIBRARIES}")
set_property(TARGET ml_classification PROPERTY CXX_STANDARD 17)
