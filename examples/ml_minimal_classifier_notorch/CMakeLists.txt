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
find_package(VTK REQUIRED)
#set (Torch_DIR "/home/mlcfd/Repositories/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} ${VTK_CXX_FLAGS}")
target_link_libraries(ml_classification PRIVATE "${TORCH_LIBRARIES}")
target_link_libraries(ml_classification PRIVATE "${VTK_LIBRARIES}")
set_property(TARGET ml_classification PROPERTY CXX_STANDARD 17)