
# Add AMR executable
add_executable(cylinder_amr_generation)
set_target_properties(cylinder_amr_generation
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_link_libraries(cylinder_amr_generation PUBLIC irl)
target_link_libraries(cylinder_amr_generation PUBLIC "-lpthread")
target_link_libraries(cylinder_amr_generation PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(cylinder_amr_generation PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(cylinder_amr_generation PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(cylinder_amr_generation SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp/")
target_include_directories(cylinder_amr_generation SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/triangle/")
target_link_libraries(cylinder_amr_generation INTERFACE PUBLIC Eigen3::Eigen)
if(${IRL_USE_CGAL})
  target_link_libraries(cylinder_amr_generation INTERFACE PUBLIC CGAL)
endif()
target_include_directories(cylinder_amr_generation SYSTEM PUBLIC "")

# Define the executable target
target_sources(cylinder_amr_generation PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/amr_generation.cpp)
target_sources(cylinder_amr_generation PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/geometry.h)

# Add cylinder_verification executable
add_executable(cylinder_verification)
set_target_properties(cylinder_verification
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_include_directories(cylinder_verification PRIVATE "./")
target_link_libraries(cylinder_verification PUBLIC irl)
target_link_libraries(cylinder_verification PUBLIC "-lpthread")
target_link_libraries(cylinder_verification PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(cylinder_verification PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(cylinder_verification PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(cylinder_verification SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp/")
target_include_directories(cylinder_verification SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/triangle/")
target_link_libraries(cylinder_verification INTERFACE PUBLIC Eigen3::Eigen)
target_link_libraries(cylinder_verification INTERFACE PUBLIC CGAL)
target_include_directories(cylinder_verification SYSTEM PUBLIC "")

target_sources(cylinder_verification PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/cylinder_confirmation.cpp)
target_sources(cylinder_verification PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/geometry.h)

# Add cylinder_test_configuration executable
add_executable(cylinder_test_configuration)
set_target_properties(cylinder_test_configuration
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_link_libraries(cylinder_test_configuration PUBLIC irl)
target_link_libraries(cylinder_test_configuration PUBLIC "-lpthread")
target_link_libraries(cylinder_test_configuration PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(cylinder_test_configuration PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(cylinder_test_configuration PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(cylinder_test_configuration SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp/")
target_include_directories(cylinder_test_configuration SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/triangle/")
target_link_libraries(cylinder_test_configuration INTERFACE PUBLIC Eigen3::Eigen)
target_link_libraries(cylinder_test_configuration INTERFACE PUBLIC CGAL)
target_include_directories(cylinder_test_configuration SYSTEM PUBLIC "")

target_sources(cylinder_test_configuration PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/test_config.cpp)
target_sources(cylinder_test_configuration PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/geometry.h)

# Add timing executable
add_executable(cylinder_timing)
set_target_properties(cylinder_timing
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_link_libraries(cylinder_timing PUBLIC irl)
target_link_libraries(cylinder_timing PUBLIC "-lpthread")
target_link_libraries(cylinder_timing PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(cylinder_timing PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(cylinder_timing PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(cylinder_timing SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp/")
target_include_directories(cylinder_timing SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/triangle/")
target_link_libraries(cylinder_timing INTERFACE PUBLIC Eigen3::Eigen)
target_link_libraries(cylinder_timing INTERFACE PUBLIC CGAL)
target_include_directories(cylinder_timing SYSTEM PUBLIC "")

target_sources(cylinder_timing PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/irl_timing.cpp)
target_sources(cylinder_timing PRIVATE ${IRL_TOOLS_SRC_DIR}/cylinder_testing/geometry.h)

add_custom_target(cylinder_testing
    DEPENDS cylinder_amr_generation cylinder_verification cylinder_test_configuration cylinder_timing
)
