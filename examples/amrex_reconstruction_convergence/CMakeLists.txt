set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  amrex_reconstruction_convergence
  HDRS
  SRCS
  "main.cpp"
  "../implicit_surface_reconstruction/basic_mesh.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/amrex_reconstruction_convergence"
  COLLECTION
  "examples_amrex_reconstruction_convergence"
)
