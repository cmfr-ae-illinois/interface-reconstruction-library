IRLAddExample_cpp(
  TARGET_NAME
  initialize_vf
  HDRS
  "functions.h"
  "surfaces.h"
  SRCS
  "main.cpp"
  "functions.cpp"
  "../variant_advector/basic_mesh.cpp"
  "output.cpp"

  # implicit surface definitions
  "./surface_definitions/sphere.cpp"
  "./surface_definitions/ellipsoid.cpp"
  "./surface_definitions/torus.cpp"
  "./surface_definitions/genus.cpp"
  "./surface_definitions/orthocircle.cpp"
  "./surface_definitions/wineglass.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/initialize_vf"
  COLLECTION
  "examples_initialize_vf"
)