# sandbox/CMakeLists.txt

# Define the executable target
add_executable(litle_test litle_test.cpp)

# Link the libraries from the project (adjust as needed)
target_link_libraries(litle_test PRIVATE irl Eigen3::Eigen CGAL)

# Include directories from the project
target_include_directories(litle_test PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/external/abseil-cpp)
