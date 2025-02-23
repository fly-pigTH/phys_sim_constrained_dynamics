file(GLOB_RECURSE link_cpp_files "./src/*.cpp")

add_library(link SHARED ${link_cpp_files})
target_link_libraries(link basic)

# Add a test project.
add_executable(test_link "./main.cpp")
target_link_libraries(test_link link)