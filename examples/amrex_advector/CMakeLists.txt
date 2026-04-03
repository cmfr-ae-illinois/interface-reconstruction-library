set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  amrex_advector
  HDRS
  "vtk.h"
  "amrcore_advector.h"
  "amrcore_advector.tpp"
  "deformation_3d.h"
  SRCS
  "main.cpp"
  "vtk.cpp"
  # "amrcore_advector.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/amrex_advector"
  COLLECTION
  "examples_amrex_advector"
)