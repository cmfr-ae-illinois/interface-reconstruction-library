set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  sphere_initialization
  HDRS
  "basic_mesh.h"
  "data.h"
  "reconstruction_types.h"
  "solver.h"
  "vtk.h"
  SRCS
  "main.cpp"
  "basic_mesh.cpp"
  "solver.cpp"
  "reconstruction_types.cpp"
  "vtk.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/sphere_initialization"
  COLLECTION
  "examples_sphere_initialization"
)