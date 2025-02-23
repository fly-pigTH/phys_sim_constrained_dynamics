file(GLOB_RECURSE contact_cpp_files "./src/*.cpp")

add_library(contact SHARED ${contact_cpp_files})
target_link_libraries(contact basic)

# Add a test project.
add_executable(test_contact "./main.cpp")
target_link_libraries(test_contact contact)