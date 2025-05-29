IRLAddExample_cpp(
  TARGET_NAME
    marching_squares
  HDRS
    "irl2d.h"
    "data.h"
    "basic_mesh.h"
  SRCS
    "main.cpp"
    "irl2d.cpp"
    "basic_mesh.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/marching_squares"
  COLLECTION
    "examples_marching_squares"
)

