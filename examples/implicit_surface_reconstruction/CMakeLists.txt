set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  implicit_surface_reconstruction
  HDRS
  "basic_mesh.h"
  "data.h"
  "reconstruction_types.h"
  "solver.h"
  "vtk.h"
  "initialization.h"
  "binary.h"
  "reconstruction.h"
  "surface_select.h"
  SRCS
  "main.cpp"
  "basic_mesh.cpp"
  "solver.cpp"
  "reconstruction_types.cpp"
  "vtk.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/implicit_surface_reconstruction"
  COLLECTION
  "examples_implicit_surface_reconstruction"
)