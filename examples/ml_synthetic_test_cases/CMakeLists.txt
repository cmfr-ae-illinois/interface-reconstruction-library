IRLAddExample_cpp(
  TARGET_NAME
    ml_data_vtk
  HDRS
  SRCS
    "main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/ml_synthetic_test_cases"
  COLLECTION
    "examples_ml_synthetic_test_cases"
)
find_package(VTK REQUIRED)
#set (Torch_DIR "/home/mlcfd/Repositories/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} ${VTK_CXX_FLAGS}")
target_link_libraries(ml_data_vtk PRIVATE "${TORCH_LIBRARIES}")
target_link_libraries(ml_data_vtk PRIVATE "${VTK_LIBRARIES}")
set_property(TARGET ml_data_vtk PROPERTY CXX_STANDARD 17)