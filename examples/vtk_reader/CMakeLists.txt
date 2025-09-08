IRLAddExample_cpp(
  TARGET_NAME
    vtk_reader
  HDRS
  SRCS
    "main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/vtk_reader"
  COLLECTION
    "examples_vtk_reader"
)

find_package(VTK REQUIRED)
# include(${VTK_USE_FILE})
target_link_libraries(vtk_reader PRIVATE "${VTK_LIBRARIES}")
set (Torch_DIR "/home/mlcfd/Repositories/libtorch/share/cmake/Torch")
find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
target_link_libraries(vtk_reader PRIVATE "${TORCH_LIBRARIES}")
set_property(TARGET vtk_reader PROPERTY CXX_STANDARD 17)
