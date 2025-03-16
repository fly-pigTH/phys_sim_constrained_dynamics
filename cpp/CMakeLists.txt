cmake_minimum_required(VERSION 3.20)

project(phys_sim_constrained_dynamics)

# Project-wide setup.
set(CMAKE_BUILD_TYPE Release)   # 构建类型为release
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)      # 设置cpp标准为cpp17

# Build external libraries.  添加外部库 polyscope 的子目录，并指定其构建目录
add_subdirectory("${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/external/polyscope"
    "${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/build/polyscope")

# Configure include and lib paths. 配置项目的包含路径，使得编译器能够找到头文件
include_directories("${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/cpp")
include_directories("${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/external")
include_directories("${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/external/eigen")
include_directories("${PHYS_SIM_CONSTRAINED_DYNAMICS_ROOT}/external/polyscope/include")

# Compile each library. 添加子目录并构建子模块。注释部分描述了各个模块的依赖关系
# Dependency map:
# basic <- empty;
add_subdirectory("basic")
# link <- basic;
add_subdirectory("link")
# joint <- link, basic;
add_subdirectory("joint")
# contact <- link, basic;
add_subdirectory("contact")
# sim <- contact, joint, link, basic;
add_subdirectory("sim")

# Compile the main program.
add_executable(main "main.cpp")
target_link_libraries(main sim polyscope)   # 链接 main 可执行文件所需的库，包括 sim 和 polyscope