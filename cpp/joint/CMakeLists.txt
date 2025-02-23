file(GLOB_RECURSE joint_cpp_files "./src/*.cpp")

add_library(joint SHARED ${joint_cpp_files})
target_link_libraries(joint link basic)

# Add a test project.
add_executable(test_joint "./main.cpp")
target_link_libraries(test_joint joint)