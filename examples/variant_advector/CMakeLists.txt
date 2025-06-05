set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  variant_advector
  HDRS
  "basic_mesh.h"
  "data.h"
  "translation_3d.h"
  "reconstruction_types.h"
  "solver.h"
  "vof_advection.h"
  "vtk.h"
  SRCS
  "main.cpp"
  "basic_mesh.cpp"
  "vof_advection.cpp"
  "solver.cpp"
  "reconstruction_types.cpp"
  "translation_3d.cpp"
  "vtk.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/variant_advector"
  COLLECTION
  "examples_variant_advector"
)
