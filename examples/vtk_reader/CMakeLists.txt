find_package(VTK QUIET)
find_package(Torch QUIET)
if(VTK_FOUND AND Torch_FOUND)
message(STATUS "VTK and Torch were found; enabling support for vtk_reader example.")
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
target_link_libraries(vtk_reader PRIVATE "${VTK_LIBRARIES}")
target_link_libraries(vtk_reader PRIVATE "${TORCH_LIBRARIES}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
set_property(TARGET vtk_reader PROPERTY CXX_STANDARD 17)
else()
message(STATUS "VTK and/or Torch were not found; disabling support for vtk_reader example.")
endif()