find_package(MPI REQUIRED COMPONENTS CXX)

IRLAddExample_cpp(
  TARGET_NAME level_set_reconstruction
  HDRS 
    "level_set.h" 
    "pu_field.h" 
    "curvature.h" 
    "pu_ppic.h" 
    "visualization.h" 
    "weights/wu2.h" 
    "weights/wu4.h" 
    "weights/wendland6.h" 
    "weights/wendland4.h" 
    "diagnostics_writer.h"
  SRCS
    "main.cpp"
    "visualization.tpp"
    "../variant_advector/basic_mesh.cpp"
    "../variant_advector/reconstruction_types.cpp"
    "../variant_advector/solver.cpp"
    "../variant_advector/vof_advection.cpp"
    "../variant_advector/vtk.cpp"
    "weights/wu2.tpp"
    "weights/wu4.tpp"
    "weights/wendland6.tpp"
    "weights/wendland4.tpp"
  OUTPUT_DIR "${IRL_EXAMPLES_OUTPUT_DIR}/level_set_reconstruction"
  COLLECTION examples_level_set_reconstruction
)
target_link_libraries(level_set_reconstruction PRIVATE MPI::MPI_CXX)

if(BUILD_TESTING)
  add_executable(level_set_pu_test pu_field_test.cpp ../variant_advector/basic_mesh.cpp)
  target_link_libraries(level_set_pu_test PRIVATE irl)
  add_test(NAME level_set_pu_test COMMAND level_set_pu_test)
endif()
