IRLAddExample_cpp(
  TARGET_NAME
    spline_example
  HDRS
    "Spline.h"
  SRCS
    "main.cpp"
    "Spline.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/spline_example"
  COLLECTION
    "examples_spline"
)

