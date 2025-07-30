#including json
# Download and unpack nljson at configure time
message(STATUS "Downloading json package")
if(EXISTS "${CMAKE_SOURCE_DIR}/external/json")
    message(STATUS "json directory exists. Reseting")
    execute_process(COMMAND "git" "reset" "--hard" "-q"
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external/json)
else()
    execute_process(COMMAND "git" "clone" "https://github.com/nlohmann/json.git"
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external/)
endif()


if(result)
    message(FATAL_ERROR "CMake step for nljson failed: ${result}")
endif()
# Add nljson directly to our build. This defines
# the gtest and gtest_main targets.


add_executable(paraboloid_amr_generate_result)
set_target_properties(paraboloid_amr_generate_result
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_include_directories(paraboloid_amr_generate_result PRIVATE "./")
target_link_libraries(paraboloid_amr_generate_result PUBLIC irl)
target_link_libraries(paraboloid_amr_generate_result PUBLIC "-lpthread")
target_link_libraries(paraboloid_amr_generate_result PUBLIC MPI::MPI_CXX)
target_link_libraries(paraboloid_amr_generate_result PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(paraboloid_amr_generate_result PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(paraboloid_amr_generate_result PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(paraboloid_amr_generate_result SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/abseil-cpp/")
target_include_directories(paraboloid_amr_generate_result SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/triangle/")
target_include_directories(paraboloid_amr_generate_result SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/json/include")
target_link_libraries(paraboloid_amr_generate_result INTERFACE PUBLIC Eigen3::Eigen)
if(${IRL_USE_CGAL})
  target_link_libraries(paraboloid_amr_generate_result INTERFACE PUBLIC CGAL)
endif()
target_include_directories(paraboloid_amr_generate_result SYSTEM PUBLIC "")
target_include_directories(paraboloid_amr_generate_result SYSTEM PUBLIC "${MPI_INCLUDE_PATH}")

target_sources(paraboloid_amr_generate_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/amr_main.cpp")
target_sources(paraboloid_amr_generate_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/amr_run.cpp")
target_sources(paraboloid_amr_generate_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/amr_run.h")

# Add IRL executable
add_executable(paraboloid_irl_confirm_result)
set_target_properties(paraboloid_irl_confirm_result
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMPILE_FLAGS "${CMAKE_CXX_FLAGS}"
  INTERPROCEDURAL_OPTIMIZATION False
  CXX_STANDARD 17)

target_include_directories(paraboloid_irl_confirm_result PRIVATE "./")
target_link_libraries(paraboloid_irl_confirm_result PUBLIC irl)
target_link_libraries(paraboloid_irl_confirm_result PUBLIC "-lpthread")
target_link_libraries(paraboloid_irl_confirm_result PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(paraboloid_irl_confirm_result PUBLIC "-lquadmath")

if(APPLE)
  target_link_libraries(paraboloid_irl_confirm_result PUBLIC "-framework CoreFoundation")
endif()

target_include_directories(paraboloid_irl_confirm_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/abseil-cpp/")
target_include_directories(paraboloid_irl_confirm_result SYSTEM PUBLIC "${IRL_ROOT_LOCATION}/external/triangle/")
target_include_directories(paraboloid_irl_confirm_result SYSTEM PUBLIC "${PROJECT_SOURCE_DIR}/external/json/include")
target_link_libraries(paraboloid_irl_confirm_result INTERFACE PUBLIC Eigen3::Eigen)
target_link_libraries(paraboloid_irl_confirm_result INTERFACE PUBLIC CGAL)
target_include_directories(paraboloid_irl_confirm_result SYSTEM PUBLIC "")


target_sources(paraboloid_irl_confirm_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/irl_main.cpp")
target_sources(paraboloid_irl_confirm_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/irl_run.h")
target_sources(paraboloid_irl_confirm_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/irl_run.tpp")
target_sources(paraboloid_irl_confirm_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/my_dodecahedron.h")
target_sources(paraboloid_irl_confirm_result PRIVATE "${IRL_TOOLS_SRC_DIR}/paraboloid_testing/my_dodecahedron.tpp")

add_custom_target(paraboloid_testing
    DEPENDS paraboloid_amr_generate_result paraboloid_irl_confirm_result
)
