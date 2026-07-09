IRLAddExample_fortran(
  TARGET_NAME
    f_cutting_methods
  SRCS
    "cutting_methods.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

IRLAddExample_fortran(
  TARGET_NAME
    f_localized_separator_link
  SRCS
    "localized_separator_link.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

IRLAddExample_fortran(
  TARGET_NAME
    f_reconstruction_elvira
  SRCS
    "reconstruction_elvira.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

IRLAddExample_fortran(
  TARGET_NAME
    f_reconstruction_mof
  SRCS
    "reconstruction_mof.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

IRLAddExample_fortran(
  TARGET_NAME
    f_cube_cutting
  SRCS
    "cube_cutting.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

IRLAddExample_fortran(
  TARGET_NAME
    f_polygon_cutting
  SRCS
    "polygon_cutting.f90"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)

find_package(VTK QUIET)
if(VTK_FOUND)
message(STATUS "VTK was found; enabling VTK support for Fortran ml_classifier examples.")
IRLAddExample_cpp(
  TARGET_NAME
    f_ml_classifier
  HDRS
    "ml_classifier/ml_classifier_c_api.h"
  SRCS
    "ml_classifier/precision_dp.f90"
    "ml_classifier/ml_classifier.f90"
    "ml_classifier/ml_classifier_c_api.f90"
    "ml_classifier/main.cpp"
  OUTPUT_DIR
    "${IRL_EXAMPLES_OUTPUT_DIR}/fortran"
  COLLECTION
    "examples_fortran"
)
target_link_libraries(f_ml_classifier PRIVATE "${VTK_LIBRARIES}")
else()
message(STATUS "VTK was not found; disabling VTK support for Fortran ml_classifier examples.")
endif()
